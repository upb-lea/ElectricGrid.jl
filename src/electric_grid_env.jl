mutable struct ElectricGridEnv <: AbstractEnv
    verbosity
    nc
    sys_d
    action_space
    state_space
    done
    inner_featurize
    featurize
    prepare_action
    reward_function
    x0
    x
    t0
    t
    ts
    state
    maxsteps
    steps
    state_ids
    v_dc
    v_dc_arr
    norm_array
    convert_state_to_cpu
    reward
    action
    action_ids
    action_delay_buffer
    A
    B
    C
    D
    state_parameters
    y #holds all of the inductor voltages and capacitor currents
    agent_dict
end


"""
    ElectricGridEnv(...)

# Description
Returns an environment consisting of an electrical power grid as control plant,
which can interact via this interface with a reinforcement learing agent from
https://juliareinforcementlearning.org/

# Arguments

# Keyword Arguments
- `maxsteps::Int`: The number of time steps that the simulation is run.
- `ts::Float`: Sampling time by which the environment is evolved per step.
- `action_space::Space`: Defines the valide space per action.
- `state_space::Space`: Defines the valide space per state. Default is [-1, 1] since the
    states of the env will be normalized.
- `prepare_action::function(env::ElectricGridEnv)`: Function to adjust, change, prepare the actions
    before they are applied to the env during a step(). Per default it returns the action
    without changes.
- `featurize::function(x:Vector, t:Float; env::ElectricGridEnv, name::String)`: Function to adjust
   the state before it is returned by the env. For example here
   reference values can be added to provide the to an RL agent or other feature engineering
   to improve the learning.
- `reward_function::function(env::ElectricGridEnv, name::String) )`: Function to define the reward
    for a (named) policy. Return 0 per default.
- `CM::Matrix`: Conectivity matrix to define the structure of the electric power grid
    (for more details see NodeConstructor).
- `num_sources::Int`: Number of sources in the electric power grid.
- `num_loads::Int`: Number of loads in the electric power grid.
- `parameter::Dict`: Dictonary to define the parameter of the grid. Entries can be "grid",
    "source", "load", "cable". Here, e.g. the electric components are defined.
- `x0::Vector`: Initial states which will be used in reset!().
- `t0::Float`: Initial time where the env starts from simulating.
- `state_ids::Vector{String}`: IDs of the states. Can be used for plotting or collecting
    trajectories during a simulation with the DataHook.
- `convert_state_to_cpu::Bool`: Flag if the states are converted to cpu if calculated one
    gpu before
- `use_gpu::Bool`: Flag if the simulation is done on gpu (if possible).
- `reward::Float`: Reward which the action achieved, applied in the env.
- `action::Vector{Float}`: Action which is applied to the env.
- `action_ids::Vector{String}`: IDs of the actions. Can be used for plotting or collecting
    trajectories during a simulation with the DataHook.
- `state_action_delay::Int`: Specifies the number of simulation steps it is delayed from
    giving an action to the env to when it is applied to the plant.
- `t_end::Float`: Spezifies the end time of an episode. If defined it overwrites maxsteps.
- `verbosity::Int`: Defines the level of how much information (warnings, info,...) is
    printed to the console (0 - nothing, 1 - warn, 2 - debug(warn, info,...)).
- `state_ids_RL::Vector{String}`: State IDs which are given to RL agents.
- `action_ids_RL::Vector{String}`: Action IDs which are given to RL agents.

# Return Values
- `ElectricGridEnv::ElectricGridEnv`: Environment an agent can interact with


"""
function ElectricGridEnv(;
    maxsteps=500,
    ts=1 / 10_000,
    action_space=nothing,
    state_space=nothing,
    prepare_action=nothing,
    featurize=nothing,
    reward_function=nothing,
    CM=nothing,
    num_sources=nothing,
    num_loads=nothing,
    parameters=nothing,
    x0=nothing,
    t0=0.0,
    state_ids=nothing,
    convert_state_to_cpu=true,
    use_gpu=false,
    reward=nothing,
    action=nothing,
    action_ids=nothing,
    action_delay=1,
    t_end=nothing,
    verbosity=0,
    agent_dict=nothing
)

    if !(isnothing(t_end))
        maxsteps = floor(t_end / ts) + 1
    end

    if !(isnothing(num_sources) || isnothing(num_loads))
        nc = NodeConstructor(num_sources=num_sources, num_loads=num_loads, CM=CM,
            parameters=parameters, ts=ts, verbosity=verbosity)
    else

        if isnothing(parameters)
            # Construct standard env with 2 sources, 1 load
            @info("Three phase electric power grid with 2 sources and 1 load is created!
            Parameters are drawn randomly! To change, please define parameters
            (see node_constructor)")
            CM = [0.0 0.0 1.0
                0.0 0.0 2
                -1.0 -2.0 0.0]
            nc = NodeConstructor(num_sources=2, num_loads=1, CM=CM, verbosity=verbosity)
        else
            if haskey(parameters, "source") && haskey(parameters, "load")
                num_sources = length(parameters["source"])
                num_loads = length(parameters["load"])
            elseif haskey(parameters, "source")
                num_sources = length(parameters["source"])
                num_loads = 0
            else
                @info("Three phase electric power grid with 2 sources and 1 load is created!
                 Parameters are drawn randomly! To change, please define parameters
                 (see node_constructor)")
                num_sources = 2
                num_loads = 1
            end
            nc = NodeConstructor(num_sources=num_sources, num_loads=num_loads, CM=CM,
                parameters=parameters, ts=ts, verbosity=verbosity)
        end
    end

    A, B, C, D = GetSystem(nc)
    Ad = exp(A * ts) #fastExpm(A*ts) might be a better option
    #Bd = A \ (Ad - I) * B #This may be bad for large sizes, maybe QR factorise, then use ldiv!
    Bd = (Ad - I) * B #
    ldiv!(factorize(A), Bd)
    sys_d = HeteroStateSpace(Ad, Bd, C, D, Float64(ts))
    state_parameters = GetStateParameters(nc)

    if use_gpu
        Ad = CuArray(A)
        Bd = CuArray(B)
        C = CuArray(C)
        if isa(D, Array)
            D = CuArray(D)
        end
    end

    if isnothing(action_space)
        action_space = Space([-1.0 .. 1.0 for i = 1:length(sys_d.B[1, :])],)
    end

    if isnothing(state_ids)
        if isnothing(nc)
            state_ids = []
            @warn("No state_ids array specified - observing states with DataHook not
                    possible")
        else
            state_ids = GetStateIds(nc)
        end
    end

    if isnothing(action_ids)
        if isnothing(nc)
            action_ids = []
            @warn("No state_ids array specified - observing states with DataHook not
                    possible")
        else
            action_ids = GetActionIds(nc)
        end
    end

    if isnothing(agent_dict)
        agent_dict = Dict()
        for ns in 1:nc.num_sources
            if nc.parameters["source"][ns]["control_type"] == "RL"
                if nc.parameters["source"][ns]["mode"] == "ElectricGrid_ddpg"
                    name = nc.parameters["source"][ns]["mode"] * "_$ns"
                else
                    name = nc.parameters["source"][ns]["mode"]
                end
                ssa = "source$ns"
                if haskey(agent_dict, name)
                    push!(agent_dict[name]["source_number"], ns)
                    append!(agent_dict[name]["state_ids"], filter(x -> split(x, "_")[1] == ssa, state_ids)[:])
                    append!(agent_dict[name]["action_ids"], filter(x -> split(x, "_")[1] == ssa, action_ids)[:])
                else
                    agent_dict[name] = Dict(
                        "source_number" => [ns],
                        "mode" => nc.parameters["source"][ns]["mode"],
                    )
                    agent_dict[name]["state_ids"] = filter(x -> split(x, "_")[1] == ssa, state_ids)
                    agent_dict[name]["action_ids"] = filter(x -> split(x, "_")[1] == ssa, action_ids)
                end
            end
        end
    end


    inner_featurize = function (x0=nothing, t0=nothing; env=nothing, name=nothing)
        if !isnothing(name)
            if name == "classic"
                return env.state
            else
                state = env.state[findall(x -> x in env.agent_dict[name]["state_ids"], env.state_ids)]
                if !isnothing(featurize)
                    state = featurize(state, env, name)
                end
                return state
            end
        elseif isnothing(env)
            return x0
        else
            return env.state
        end
    end


    if isnothing(prepare_action)
        prepare_action = function (env)
            env.action
        end
    end

    if isnothing(reward_function)
        reward_function = function (env, name=nothing)
            return 0.0
        end
    end

    if isnothing(x0)
        x0 = [0.0 for i = 1:length(sys_d.A[1, :])]
    end

    x = x0
    t = t0
    state = inner_featurize(x0, t0)

    if isnothing(state_space)
        state_space = Space([-1.0 .. 1.0 for i = 1:length(state)],)
    end

    if use_gpu
        if isa(x0, Array)
            x0 = CuArray(x0)
        end
        if !(convert_state_to_cpu) && isa(state, Array)
            state = CuArray(state)
        end
    end

    vdc_fixed = 0
    v_dc = ones(nc.num_sources)  # vector to store evaluated v_dc_arr (constants and
    # functions) in the env, needed e.g. in the DataHook
    v_dc_arr = []  # array to store all functions for v_dc as well as constants
    for (source_number, source) in enumerate(nc.parameters["source"])
        if haskey(source, "source_type")
            if source["source_type"] == "ideal"
                # set vdc for that source
                if haskey(source, "vdc")
                    v_dc[source_number] = source["vdc"]
                    fun = (env, G, T) -> source["vdc"]
                    push!(v_dc_arr, fun)
                else
                    v_dc[source_number] = 800
                    fun = (env, G, T) -> 800
                    push!(v_dc_arr, fun)
                    vdc_fixed += 1
                end
            elseif source["source_type"] == "pv"
                #v_dc[source_number] = :(GetV(SolarArray, env.x[1]*env.action, G, T))
                #TODO : how to calculate i_dc in 3-phase grid? Which current of env to use?
                #TODO : $source_number does only fit if all L filters! Otherwise how to
                #       define the offet for $source_number?!?!?
                # TODO built pv module from parameter dict - where to define? In env?
                ModulePV = SolarModule()
                SolarArr = SolarArray(; ModuleParameters=ModulePV)
                # find(x -> .... source$source_number_i_L in state_ids)
                fun = (env, G, T) -> GetV(:($SolarArr),
                    env.x[:($source_number)] * env.action, G, T)
                push!(v_dc_arr, fun)

                # first value set to 0
                v_dc[source_number] = 0
            else
                @warn "sourceType not known! vdc set to fixed value"
                v_dc[source_number] = 800
                fun = (env, G, T) -> 800
                push!(v_dc_arr, fun)
                vdc_fixed += 1
            end
        else
            @warn "sourceType not defined! vdc set to fixed value, if not wanted please
                   define nc.parameters -> source -> source_type (e.g. = ideal)"
            v_dc[source_number] = 800
            fun = (env, G, T) -> 800
            push!(v_dc_arr, fun)
            vdc_fixed += 1
        end
    end
    vdc_fixed > 0 && @warn("$vdc_fixed DC-link voltages set to 800 V - please define in
                            nc.parameters -> source -> vdc.")

    if verbosity > 1
        @info "Normalization is done based on the defined parameter limits."
        @info "Time simulation run time: $((maxsteps - 1)*ts) [s] ~> $(Int(maxsteps)) steps"
    end
    states = GetStateIds(nc)

    i_limit_fixed = 0
    v_limit_fixed = 0
    norm_array = ones(length(states))

    for (source_number, source) in enumerate(nc.parameters["source"])
        # set norm_array based on in parameters defined limits
        for state_index in GetSourceStateIndices(nc, [source_number])["source$source_number"]["state_indices"]
            if contains(states[state_index], "_i")
                if haskey(source, "i_limit")
                    norm_array[state_index] = source["i_limit"]
                else
                    i_limit_fixed += 1
                    norm_array[state_index] = 1.15 * sqrt(2) *
                                              nc.parameters["source"][source_number]["pwr"] /
                                              (3 * nc.parameters["grid"]["v_rms"])
                    nc.parameters["source"][source_number]["i_limit"] = norm_array[state_index]
                end
            elseif contains(states[state_index], "_v")
                if haskey(source, "v_limit")
                    norm_array[state_index] = source["v_limit"]
                else
                    v_limit_fixed += 1
                    norm_array[state_index] = 1.5 * nc.parameters["source"][source_number]["vdc"]
                    nc.parameters["source"][source_number]["v_limit"] = norm_array[state_index]
                end
            end
        end
    end

    for (load_number, load) in enumerate(nc.parameters["load"])
        for state_index in GetLoadStateIndices(nc, [load_number])["load$load_number"]["state_indices"]
            if contains(states[state_index], "_i")
                if haskey(load, "i_limit")
                    norm_array[state_index] = load["i_limit"]
                else
                    i_limit_fixed += 1
                    norm_array[state_index] = 1000.0
                end
            elseif contains(states[state_index], "_v")
                if haskey(load, "v_limit")
                    norm_array[state_index] = load["v_limit"]
                else
                    v_limit_fixed += 1
                    norm_array[state_index] = 1.15 * nc.parameters["grid"]["v_rms"] * sqrt(2)
                end
            end
        end
    end

    for (cable_number, cable) in enumerate(nc.parameters["cable"])
        for state_index in GetCableStateIndices(nc, [cable_number])["cable$cable_number"]["state_indices"]
            if contains(states[state_index], "_i")
                if haskey(cable, "i_limit")
                    norm_array[state_index] = cable["i_limit"]
                else
                    i_limit_fixed += 1
                    norm_array[state_index] = 1000
                end
            elseif contains(states[state_index], "_v")
                if haskey(cable, "v_limit")
                    norm_array[state_index] = cable["v_limit"]
                else
                    v_limit_fixed += 1
                    norm_array[state_index] = 1.15 * nc.parameters["grid"]["v_rms"] * sqrt(2)
                end
            end
        end
    end

    if isnothing(reward)
        reward = 0.0
    end

    if isnothing(action)
        action = zeros(length(action_space))
    end

    if action_delay == 0
        action_delay_buffer = nothing
    else
        action_delay_buffer = CircularBuffer{Vector{Float64}}(action_delay)
        fill!(action_delay_buffer, zeros(length(action_space)))
    end

    y = (A * Vector(x) + B * (Vector(action))) .* (state_parameters)

    ElectricGridEnv(verbosity, nc, sys_d, action_space, state_space,
        false, inner_featurize, featurize, prepare_action, reward_function,
        x0, x, t0, t, ts, state, maxsteps, 0, state_ids,
        v_dc, v_dc_arr, norm_array, convert_state_to_cpu,
        reward, action, action_ids, action_delay_buffer,
        A, B, C, D, state_parameters, y, agent_dict)
end

RLBase.action_space(env::ElectricGridEnv) = env.action_space
RLBase.state_space(env::ElectricGridEnv) = env.state_space
RLBase.reward(env::ElectricGridEnv) = env.reward
RLBase.DynamicStyle(env::ElectricGridEnv) = RLBase.Simultaneous

function RLBase.reward(env::ElectricGridEnv, name::String)
    return env.reward_function(env, name)
end

RLBase.is_terminated(env::ElectricGridEnv) = env.done
RLBase.state(env::ElectricGridEnv) = env.state

function RLBase.state(env::ElectricGridEnv, name::String)
    return env.inner_featurize(; env=env, name=name)
end

"""
    reset!(env)

Resets the environment. The state is set to x0, the time to t0, reward to zero and done to
false.
"""
function RLBase.reset!(env::ElectricGridEnv)
    env.state = env.convert_state_to_cpu ? Array(env.inner_featurize(env.x0, env.t0)) : env.inner_featurize(env.x0, env.t0)
    env.x = env.x0
    env.y = fill!(env.y, 0.0) #TODO: y0
    if !isnothing(env.action_delay_buffer)
        empty!(env.action_delay_buffer)
        fill!(env.action_delay_buffer, zeros(length(env.action_space)))
        env.action = zeros(length(env.action_space))
    end
    env.t = env.t0
    env.steps = 0
    env.reward = 0.0
    env.done = false
    return nothing
end

action_space(env::ElectricGridEnv) = RLBase.action_space(env)
state_space(env::ElectricGridEnv) = RLBase.state_space(env)
reward(env::ElectricGridEnv) = RLBase.reward(env)
DynamicStyle(env::ElectricGridEnv) = RLBase.DynamicStyle(env)
reward(env::ElectricGridEnv, name::String) = RLBase.reward(env, name)
is_terminated(env::ElectricGridEnv) = RLBase.is_terminated(env)
state(env::ElectricGridEnv) = RLBase.state(env)
state(env::ElectricGridEnv, name::String) = RLBase.state(env, name)
reset!(env::ElectricGridEnv) = RLBase.reset!(env)


"""
    env(action)

# Description
Basic function to interact with an RL agent or other controller. Takes the action and evolvs
the system defined by the linear state-space system using that action.

# Keyword Arguments
- `action::Vector`: Action which is applied to the env. One per source per phase. Will be
    multiplied with vdc and can be interpreted as modulation index.

"""
function (env::ElectricGridEnv)(action)
    env.steps += 1
    t = [env.t, env.t + env.ts]
    env.t = t[2]

    if !isnothing(env.action_delay_buffer)
        env.action = env.action_delay_buffer[1]
        push!(env.action_delay_buffer, action)
    else
        env.action = action
    end

    # TODO V2: define G and T via data_set or stochastic process next to SolarArray
    G = 1000
    T = 27
    env.v_dc = [vdc(env, G, T) for vdc in env.v_dc_arr]
    env.action = env.action .* repeat(env.v_dc / 2, inner=env.nc.parameters["grid"]["phase"])

    env.action = env.prepare_action(env)

    if env.sys_d.A isa CuArray && !(env.action isa CuArray)
        if env.action isa Array
            env.action = CuArray(env.action)
        else
            env.action = CuArray([env.action])
        end
    end
    u = [env.action env.action]

    xout_d = CustomLsim(env.sys_d, u, t, x0=env.x)
    env.x = xout_d[:, 2]

    if env.convert_state_to_cpu
        env.state = Array(xout_d)'[2, :] ./ env.norm_array
    else
        env.state = xout_d'[2, :] ./ env.norm_array
    end

    env.state = env.inner_featurize(; env=env)
    env.reward = env.reward_function(env)
    env.done = (env.steps >= env.maxsteps) || (any(abs.(env.x ./ env.norm_array) .> 1))

    if env.done
        if any(abs.(env.x ./ env.norm_array) .> 1)
            states_exceeded = findall(abs.(env.x ./ env.norm_array) .> 1)
            if env.verbosity > 0
                @warn ("The state(s) $(env.state_ids[states_exceeded]) exceeded limit(s)
                    -> episode abort")
                @warn ("Corresponding limit(s): $(env.norm_array[states_exceeded]),
                    corresponding index: $(states_exceeded)")
            end
        end
    end

    # calcultaing the inductor voltages and capacitor currents
    env.y = (env.A * Vector(env.x) + env.B * (Vector(env.action))) .* (env.state_parameters)
end

function GetVDC_PV(I)

    V_dc = I * N_cell * P_cell

end
