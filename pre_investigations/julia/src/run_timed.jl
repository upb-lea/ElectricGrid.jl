using ReinforcementLearning
using TimerOutputs

(agent::Agent)(env, timer::TimerOutput) = agent.policy(env, timer::TimerOutput)

function (p::DDPGPolicy)(env, timer::TimerOutput, player::Any = nothing)
    p.update_step += 1

    if p.update_step <= p.start_steps
        p.start_policy(env)
    else
        D = device(p.behavior_actor)
        s = DynamicStyle(env) == SEQUENTIAL ? state(env) : state(env, player)
        s = Flux.unsqueeze(s, ndims(s) + 1)
        s = send_to_device(D, s)

        actions = p.behavior_actor(s) |> vec

        actions = actions |> send_to_host
        # TODO: use alternative to clamp - to run solely on GPU
        c = clamp.(actions .+ randn(p.rng, p.na) .* repeat([p.act_noise], p.na), -p.act_limit, p.act_limit)
        p.na == 1 && return c[1]
        c
    end
end

function (agent::Agent)(stage::AbstractStage, env::AbstractEnv, timer::TimerOutput)
    update!(agent.trajectory, agent.policy, env, stage)
    update!(agent.policy, agent.trajectory, env, stage)
end

function (agent::Agent)(stage::PreExperimentStage, env::AbstractEnv, timer::TimerOutput)
    update!(agent.policy, agent.trajectory, env, stage)
end

function (agent::Agent)(stage::PreActStage, env::AbstractEnv, action, timer::TimerOutput)
    update!(agent.trajectory, agent.policy, env, stage, action)
    update!(agent.policy, agent.trajectory, env, stage)
end

function run(policy::AbstractPolicy,
    env::AbstractEnv,
    stop_condition = StopAfterEpisode(1),
    hook = EmptyHook(),
    timer::TimerOutput)

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

            env(action)

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