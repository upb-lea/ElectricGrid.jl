using DrWatson
@quickactivate "MicroGridSimWithRL"

using GLMakie
using GraphMakie
using Graphs

g = wheel_graph(10)
f, ax, p = graphplot(g,
        edge_width = [3.0 for i in 1:ne(g)],
        edge_color = [:black for i in 1:ne(g)],
        node_size = [20 for i in 1:nv(g)],
        node_color = [:red for i in 1:nv(g)])
        
#hidedecorations!(ax); hidespines!(ax)
#ax.aspect = DataAspect()
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

f