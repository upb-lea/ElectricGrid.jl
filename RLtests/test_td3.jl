using ReinforcementLearning
using IntervalSets
# using ElectricGrid

env = PendulumEnv()

# inspect
@show env.action_space
@show env.observation_space
@show env

@show state_space(env)
@show action_space(env)

high = env.action_space.right
low = env.action_space.left
high_space = env.observation_space[3]


state_wrapped_env = StateTransformedEnv(
    env;
    state_mapping = x -> x ./ [1.0, 1.0, high_space.right],
    state_space_mapping = _ -> IntervalSets.ClosedInterval[-1.0..1.0, -1.0..1.0, -1.0..1.0],
    )
# wrap action space and force it between -1 and 1
wrapped_env = ActionTransformedEnv(
    state_wrapped_env; 
    action_mapping = x -> x / high,
    action_space_mapping = _ -> IntervalSets.ClosedInterval[-1.0..1.0],
    )

env.reward
# inspect
@show wrapped_env(0)

@show env.reward
@show action_space(wrapped_env)
@show legal_action_space(wrapped_env)
# @show act!(wrapped_env, 0)

RLBase.test_runnable!(wrapped_env)
RLBase.test_runnable!(state_wrapped_env)

RLBase.test_interfaces!(wrapped_env)
RLBase.test_interfaces!(state_wrapped_env)

 run(RandomPolicy(),
    wrapped_env,
    StopAfterStep(10),
    )

    
