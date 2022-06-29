import Base: run

using ReinforcementLearning
using TimerOutputs
using CUDA

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
        update!(agent.trajectory, agent.policy, env, stage, action)
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

function run(policy::AbstractPolicy,
    env::AbstractEnv,
    timer::TimerOutput,
    stop_condition = StopAfterEpisode(1),
    hook = EmptyHook(),
    )

    hook(PRE_EXPERIMENT_STAGE, policy, env)
    policy(PRE_EXPERIMENT_STAGE, env, timer)
    is_stop = false
    while !is_stop
        reset!(env)
        policy(PRE_EPISODE_STAGE, env, timer)
        hook(PRE_EPISODE_STAGE, policy, env)

        while !is_terminated(env) # one episode
            action = policy(env, timer)

            policy(PRE_ACT_STAGE, env, action, timer)
            hook(PRE_ACT_STAGE, policy, env, action)

            @timeit timer "Env calculation" begin
                if env.Ad isa CuArray
                    env(CuArray(action))
                else
                    env(action)
                end
            end

            policy(POST_ACT_STAGE, env, timer)
            hook(POST_ACT_STAGE, policy, env)

            if stop_condition(policy, env)
                is_stop = true
                break
            end
        end # end of an episode

        if is_terminated(env)
            policy(POST_EPISODE_STAGE, env, timer)  # let the policy see the last observation
            hook(POST_EPISODE_STAGE, policy, env)
        end
    end
    hook(POST_EXPERIMENT_STAGE, policy, env)
    hook
end