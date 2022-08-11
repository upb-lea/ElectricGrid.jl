import Base: run

using ReinforcementLearning
using CUDA
using TimerOutputs

(agent::Agent)(env, timer::TimerOutput) = agent.policy(env, timer::TimerOutput)

function (p::DDPGPolicy)(env, timer::TimerOutput, player::Any = nothing)
    p.update_step += 1

    if p.update_step <= p.start_steps
        @timeit timer "Agent calculate" begin
            if typeof(device(p)) == Val{:cpu}
                p.start_policy(env)
            else
                CUDA.@sync begin
                    p.start_policy(env)
                end
            end
        end

    else
        @timeit timer "Agent prepare data" begin
            D = device(p.behavior_actor)
            s = DynamicStyle(env) == SEQUENTIAL ? state(env) : state(env, player)
            s = Flux.unsqueeze(s, ndims(s) + 1)
            s = send_to_device(D, s)
        end
        
        @timeit timer "Agent calculate" begin
            if typeof(device(p)) == Val{:cpu}
                actions = p.behavior_actor(s) |> vec
            else
                CUDA.@sync begin
                    actions = p.behavior_actor(s) |> vec
                end
            end
        end

        @timeit timer "Agent prepare data" begin
            actions = actions |> send_to_host
            # TODO: use alternative to clamp - to run solely on GPU
            c = clamp.(actions .+ randn(p.rng, p.na) .* repeat([p.act_noise], p.na), -p.act_limit, p.act_limit)
            p.na == 1 && return c[1]
            c
        end
    end
end

function RLBase.update!(
    p::DDPGPolicy,
    traj::CircularArraySARTTrajectory,
    ::AbstractEnv,
    ::PreActStage,
    timer::TimerOutput
)
    length(traj) > p.update_after || return
    p.update_step % p.update_freq == 0 || return
    inds, batch = sample(p.rng, traj, BatchSampler{SARTS}(p.batch_size))
    update!(p, batch, timer::TimerOutput)
end

function RLBase.update!(p::DDPGPolicy, batch::NamedTuple{SARTS}, timer::TimerOutput)
    @timeit timer "policy - transfer of data to device" begin
    s, a, r, t, s′ = send_to_device(device(p), batch)
    end

    A = p.behavior_actor
    C = p.behavior_critic
    Aₜ = p.target_actor
    Cₜ = p.target_critic

    γ = p.γ
    ρ = p.ρ


    # !!! we have several assumptions here, need revisit when we have more complex environments
    # state is vector
    # action is scalar
    a′ = Aₜ(s′)
    qₜ = Cₜ(vcat(s′, a′)) |> vec
    y = r .+ γ .* (1 .- t) .* qₜ
    a = Flux.unsqueeze(a, ndims(a) + 1)
    # println("inside Policy")

    @timeit timer "gradients for Critic " begin
        gs1 = gradient(Flux.params(C)) do
            q = C(vcat(s, a)) |> vec
            loss = mean((y .- q) .^ 2)
            # ignore() do
                p.critic_loss = loss
            # end
            loss
        end
    end
    

    @timeit timer "update Critic network " begin
        update!(C, gs1)
    end
    

    @timeit timer "gradients for Actor " begin
        gs2 = gradient(Flux.params(A)) do
            loss = -mean(C(vcat(s, A(s))))
            # ignore() do
                p.actor_loss = loss
            # end
            loss
        end
    end
            

    @timeit timer "Update Actor network " begin
        update!(A, gs2)
    end


    # polyak averaging
    @timeit timer "polyak average" begin
        for (dest, src) in zip(Flux.params([Aₜ, Cₜ]), Flux.params([A, C]))
            dest .= ρ .* dest .+ (1 - ρ) .* src
        end
    end
end


function (agent::Agent)(stage::AbstractStage, env::AbstractEnv, timer::TimerOutput)
    @timeit timer "Agent trajectory" begin
        update!(agent.trajectory, agent.policy, env, stage)
    end
    @timeit timer "Agent policy update" begin
        if typeof(device(agent)) == Val{:cpu}
            update!(agent.policy, agent.trajectory, env, stage)
        else
            CUDA.@sync begin
                update!(agent.policy, agent.trajectory, env, stage)
            end
        end
    end
end

function (agent::Agent)(stage::PreExperimentStage, env::AbstractEnv, timer::TimerOutput)
    @timeit timer "Agent policy update" begin
        if typeof(device(agent)) == Val{:cpu}
            update!(agent.policy, agent.trajectory, env, stage)
        else
            CUDA.@sync begin
                update!(agent.policy, agent.trajectory, env, stage)
            end
        end
    end
end

function (agent::Agent)(stage::PreActStage, env::AbstractEnv, action, timer::TimerOutput)
    @timeit timer "Agent trajectory" begin
        # saves action in trajectory
        update!(agent.trajectory, agent.policy, env, stage, action)
    end
    @timeit timer "Agent policy update" begin
        if typeof(device(agent)) == Val{:cpu}
            update!(agent.policy, agent.trajectory, env, stage, timer)
        else
            CUDA.@sync begin
                update!(agent.policy, agent.trajectory, env, stage, timer)
            end
        end
    end
end



function run(policy::AbstractPolicy,
    env::AbstractEnv,
    timer::TimerOutput,
    stop_condition = StopAfterEpisode(1),
    hook = EmptyHook(),
    )

    @timeit timer "inside run" begin
        @timeit timer "hooks" begin
            hook(PRE_EXPERIMENT_STAGE, policy, env)
        end

        policy(PRE_EXPERIMENT_STAGE, env, timer)
        is_stop = false
        
        while !is_stop
            reset!(env)
        
            policy(PRE_EPISODE_STAGE, env, timer)
        
            @timeit timer "hooks" begin
                hook(PRE_EPISODE_STAGE, policy, env)
            end

            while !is_terminated(env) # one episode
                action = policy(env, timer)

                policy(PRE_ACT_STAGE, env, action, timer)
        
                @timeit timer "hooks" begin
                    hook(PRE_ACT_STAGE, policy, env, action)
                end

                if env.Ad isa CuArray
                    @timeit timer "Agent prepare data" begin
                        if action isa Array
                            action = CuArray(action)
                        else
                            action = CuArray([action])
                        end
                    end
                end
        
                @timeit timer "Env calculation" begin
                    env(action)
                end

                policy(POST_ACT_STAGE, env, timer)
                @timeit timer "hooks" begin
                    hook(POST_ACT_STAGE, policy, env)
                end

                if stop_condition(policy, env)
                    is_stop = true
                    break
                end
            end # end of an episode

            if is_terminated(env)
                policy(POST_EPISODE_STAGE, env, timer)  # let the policy see the last observation
                @timeit timer "hooks" begin
                    hook(POST_EPISODE_STAGE, policy, env)
                end
            end
        end
   
        @timeit timer "hooks" begin
            hook(POST_EXPERIMENT_STAGE, policy, env)
        end
        hook
    end
end