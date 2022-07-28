using CSV
using DataFrames

function get_data(episode)

    data = CSV.read("episode_data/$(episode).csv", DataFrame)
    data.state = eval.(Meta.parse.(data.state))
    data.next_state = eval.(Meta.parse.(data.next_state))
    data.action = eval.(Meta.parse.(data.action))
    
    data

end

function plotting(episode, state, NodeConstructor) # also path to data
    
    data = get_data(episode)
    
    state_list = get_states(NodeConstructor)
    
    phase = NodeConstructor.parameters["grid"]["phase"]
    fs = NodeConstructor.parameters["grid"]["fs"]
    
    t = range(0, length(data.state)-1)/fs |> collect
        
    if typeof(state) === String
        state = findall(x->x==state, state_list)
        state = state[1]
    end
    
    data_state = mapreduce(permutedims, vcat, data.state)
    
    ns = NodeConstructor.num_spp
    
    p = plot(xlabel="time", ylabel="$(state_list[state])")
    
    if phase == 1
        
        plot!(p, t, data_state[:, state+ns*0])
    
    elseif phase == 3
            
        plot!(p, t, data_state[:, state+ns*0], label="Phase 1")
        plot!(p, t, data_state[:, state+ns*1], label="Phase 2")
        plot!(p, t, data_state[:, state+ns*2], label="Phase 3")
            
    end

end