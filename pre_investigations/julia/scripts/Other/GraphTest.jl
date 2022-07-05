using DrWatson
@quickactivate "MicroGridSimWithRL"

using GLMakie
using GraphMakie
using GraphMakie.Graphs
using GraphMakie.NetworkLayout
using Makie.Colors
#using NetworkLayout

#Makie.inline!(true)

g = wheel_graph(10)
add_edge!(g, 1, 1)
f, ax, p = graphplot(g,
                        edge_width = [3.0 for i in 1:ne(g)],
                        edge_color = [colorant"black" for i in 1:ne(g)],
                        node_size = [20 for i in 1:nv(g)],
                        node_color = [colorant"red" for i in 1:nv(g)])
deregister_interaction!(ax, :rectanglezoom)

function edge_hover_action(s, idx, event, axis)
    p.edge_width[][idx]= s ? 6.0 : 3.0
    p.edge_width[] = p.edge_width[] # trigger observable
end
ehover = EdgeHoverHandler(edge_hover_action)
register_interaction!(ax, :ehover, ehover)

function node_hover_action(s, idx, event, axis)
    p.node_size[][idx] = s ? 40 : 20
    p.node_size[] = p.node_size[] # trigger observable
end
nhover = NodeHoverHandler(node_hover_action)
register_interaction!(ax, :nhover, nhover)

function node_click_action(idx, args...)
    p.node_color[][idx] = rand(RGB)
    p.node_color[] = p.node_color[]
end
nclick = NodeClickHandler(node_click_action)
register_interaction!(ax, :nclick, nclick)

function edge_click_action(idx, args...)
    p.edge_color[][idx] = rand(RGB)
    p.edge_color[] = p.edge_color[]
end
eclick = EdgeClickHandler(edge_click_action)
register_interaction!(ax, :eclick, eclick)

function node_drag_action(state, idx, event, axis)
    p[:node_pos][][idx] = event.data
    p[:node_pos][] = p[:node_pos][]
end
ndrag = NodeDragHandler(node_drag_action)
register_interaction!(ax, :ndrag, ndrag)

mutable struct EdgeDragAction
    init::Union{Nothing, Point2f} # save click position
    src::Union{Nothing, Point2f}  # save src vertex position
    dst::Union{Nothing, Point2f}  # save dst vertex position
    EdgeDragAction() = new(nothing, nothing, nothing)
end
function (action::EdgeDragAction)(state, idx, event, axis)
    edge = collect(edges(g))[idx]
    if state == true
        if action.src===action.dst===action.init===nothing
            action.init = event.data
            action.src = p[:node_pos][][edge.src]
            action.dst = p[:node_pos][][edge.dst]
        end
        offset = event.data - action.init
        p[:node_pos][][edge.src] = action.src + offset
        p[:node_pos][][edge.dst] = action.dst + offset
        p[:node_pos][] = p[:node_pos][] # trigger change
    elseif state == false
        action.src = action.dst = action.init =  nothing
    end
end
edrag = EdgeDragHandler(EdgeDragAction())
register_interaction!(ax, :edrag, edrag)

f