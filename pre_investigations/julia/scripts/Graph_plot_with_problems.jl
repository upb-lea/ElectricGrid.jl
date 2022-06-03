#using DrWatson
#@quickactivate "MicroGridSimWithRL"

using GraphMakie
using Graphs
using GLMakie
using NetworkLayout

#_______________________________________________________________________________
## Graph Operations

function Circular_Layout(p, ax)
    # gives the coordinates to the nodes such that they are laid out in a circle
    radius = 5
    # circular layout - square lattice
    for i in 1:L
        # graph coordinates
        p[:node_pos][][i] = (radius*cos((i-1)*2*π/(L)), radius*sin((i-1)*2*π/(L)))
    end
    p[:node_pos][] = p[:node_pos][]
    limits!(ax, -radius*1.5, radius*1.5, -radius*1.5, radius*1.5)
end

function AddNode(M, loc, node)

    # Adds new unconnected nodes to the system (if it is not already there)

    if size(M,1) < node
        increase = node-size(M,1)
        L_or = size(M,1)    
        L_in = L_or + increase
        new_M = Array{Float64, 2}(undef, L_in, L_in)
        new_M = fill!(new_M, 0)

        for i in 1:L_or
            for j in 1:L_or
                if i < loc && j < loc
                    new_M[i,j] = M[i,j]
                elseif i >= loc && j < loc
                    new_M[i+increase,j] = M[i,j]
                elseif i < loc && j >= loc
                    new_M[i,j+increase] = M[i,j]
                elseif i >= loc && j >= loc
                    new_M[i+increase,j+increase] = M[i,j]
                end
            end
        end

        return new_M
    else
        return M
    end
end

function AddEdge(M, node1, node2)
    # Adds a new edge to the system, and creates the nodes if they don't exist
    if size(M,1) < node1 && size(M,1) < node2
        M = AddNode(M, maximum(node1,node2))
    end

    M[node1,node2] = 1
    M[node2,node1] = 1 # bidirectional graph

    return M
end

function GetNeighbours(M, node)

    deg = sum(M[node,:])

    N = Array{Float64, 1}(undef, Int(deg))

    count = 1
    for i in 1:size(M,1)

        if M[node,i] == 1
            N[count] = i
            count = count + 1
        end
    end

    return N
end

function CircularGraph(M, L)
    #=
    In graph theory, a cycle graph Cn, sometimes simply known as an n-cycle, is
    a graph on n nodes containing a single cycle through all nodes. A different sort
    of cycle graph, termed a group cycle graph is a graph which shows cycles of a
    group as well as the connectivity between the group cycles
    =#

    M = Array{Float64, 2}(undef, L, L)
    M = fill!(M, 0)

    for i in 1:L

        if i < L
            M[i,i+1] = 1
        else
            M[1,i] = 1
        end
    end
    M = M + transpose(M)

    return M
end

function SmallWorld(M, L, Z, p)

    # Z, the short edges connecting it to its nearest neighbours, up to a distance
    # of Z/2

    # There are p x L x Z/2 shortcuts added to the network, which connect nodes
    # at random

    distance = convert(Int64, round(Z/2))
    M = Array{Float64, 2}(undef, L, L)
    M = fill!(M, 0)

    for i in 1:L

        # add the nearest neighbour connections
        # each node i should be connected to nodes mod((i - Z/2),L),... mod((i + Z/2),L)
        for j in 1:distance
            if i + j <= L
                M[i, i + j] = 1
            else
                M[i + j - L, i] = 1
            end
        end
    end

    # add the random connections
    for i in 1:convert(Int64, round(p*L*Z/2))

        node1 = convert(Int64, round((L-1)*rand() + 1))
        node2 = convert(Int64, round((L-1)*rand() + 1))

        while M[node1, node2] >= 1 || M[node2,node1] >= 1 || node1 == node2

            node1 = convert(Int64, round((L-1)*rand() + 1))
            node2 = convert(Int64, round((L-1)*rand() + 1))

            if sum(M) >= L*(L - 1)/2
                break
            end
        end

        if node2 > node1 && M[node1,node2] < 1
            M[node1,node2] = 1
        elseif node1 > node2 && M[node2,node1] < 1
            M[node2,node1] = 1
        end
    end

    M = M + transpose(M) # The electrical network may not necessarily be
    # birectional, since sources aren't capable of acting as loads.

    return M
end

function Barabasi_Albert(M, L)
    #=
    The Barabasi-Albert model is an algorithm for generating random scale-free networks
    using a preferential attachment mechanism. Several natural and human-made systems,
    including the Internet, the World-wide-web, and some social networks are thought
    to be approximately scale-free and certainly contain few nodes (called hubs) with
    unusually high degree as compared to other nodes of the network.

    The network can begin with an initial connected network of L0 nodes.

    New nodes are added to the network one at a time. Each node is connected to
    L<= L0 existing nodes with a probability that is proportional to the number of
    links that the existing nodes already have.
    =#

    if sum(M) == 0 # no connections
        if size(M,1) == 0
            M = Array{Float64, 2}(undef, L, L)
            M = fill!(M, 0)
        end
        if L > 1
            M[1,2] = 1 #just connecting the first node with the second
            M[2,1] = 1
        else
            println("The number of nodes has to be greater than 2")
        end
    end

    for i in 3:L

        for j in 1:L
            k = length(GetNeighbours(M, j))
            kt = sum(M)
            p = k/kt

            if rand() < p
                M[i,j] = 1
                M[j,i] = 1
            end
        end
    end

    return M
end

function FindpathLengthsFromNode(M, node)
    #=
    returns for each node2 in the graph the shortest distance from node to node2
    given node, returns the shortest distances from all other nodes to the given
    node. Therefore, we are returning a vector.
    =#

    #=
    Without the random long edges, the shortest path between i and j will be given
    by hopping in steps of length Z/2 along the shorter of the two arcs around the
    circle; there will be no paths of length longer than L/Z (halfway around the
    circle), and the distribution ρ(l) of path lengths will be constant for
    0 < l < L/Z. When we add shortcuts, we expect that the distribution will be
    shifted to shorter path lengths.
    =#

    #=
    We can implement a breadth-first traversal of the graph, working outward
    from node in shells. There will be a CurrentShell of nodes whose distance
    will be set to l unless they have already been visited, and nextshell which
    will be considered after the current one is finished (looking sideways
    before forward, i.e. breadth first). The idea is to sweep outward from node,
    measuring the shortest distance to every other node in the network.
    =#

    l = 0 # initialise, the distance from the node to itself is zero
    CurrentShell = [node]
    NextShell = Array{Float64, 1}(undef, 0)
    L = size(M,1)
    distances = Array{Float64, 1}(undef, L)
    distances = fill!(distances, -1.0)

    distances[node] = l

    count = 0

    while length(CurrentShell) != 0

        #=
        For each neighbour of each node in the current shell, if the distance
        to neighbour has not been set, add the node to nextShell and set the
        distance to l + 1
        =#

        for i in 1:length(CurrentShell)

            N = GetNeighbours(M, convert(Int64, round(CurrentShell[i])))

            if length(N) != 0
                for j in 1:length(N)
                    if distances[convert(Int64, round(N[j]))] == -1
                        NextShell = push!(NextShell, N[j])
                        distances[convert(Int64, round(N[j]))] = l + 1
                    end
                end
            end
        end

        CurrentShell = NextShell

        l = l + 1

        #Fail safe
        count = count + 1
        if count > L
            break
        end
    end

    return distances
end

function FindAllPathLengths(M)
    # Generates a list of all lengths (one per pair of nodes).

    L = size(M,1)
    distances = Array{Float64, 1}(undef, 0)
    distances = transpose(FindpathLengthsFromNode(M, 1)[2:L])

    for i in 2:(L-1)
        distances = [distances transpose(FindpathLengthsFromNode(M, i)[(i+1):L])]
    end

    return distances
end

#_______________________________________________________________________________
## GUI Interactions

function node_hover_action(state, idx, event, axis)
    p.node_size[][idx] = state ? 20 : 10
    p.node_size[] = p.node_size[] # trigger observable
    if p.node_color[][idx] == :red
        p.nlabels[][idx] = state ? "Source: "*repr(idx) : "" # ? is the ternary operator, closely related to the if-elseif-else
        p.nlabels[] = p.nlabels[]
        p.nlabels_distance[] = p.nlabels_distance[]
    elseif p.node_color[][idx] == :green
        p.nlabels[][idx] = state ? "Load: "*repr(idx - num_sources) : "" # ? is the ternary operator, closely related to the if-elseif-else
        p.nlabels[] = p.nlabels[]
        p.nlabels_distance[] = p.nlabels_distance[]
    end
end

function edge_hover_action(state, idx, event, axis)
    p.edge_width[][idx]= state ? 5.0 : 2.0
    p.edge_width[]  = p.edge_width[]  # trigger observable
    p.elabels[][idx] = state ? "Line: "*repr(idx) : "" # ? is the ternary operator, closely related to the if-elseif-else
    p.elabels[] = p.elabels[]
    p.elabels_distance[] = p.elabels_distance[]
end

function node_click_action(idx, args...)
    if p.node_color[][idx] == :purple
        p.node_color[][idx] = :red
    elseif p.node_color[][idx] == :red
        p.node_color[][idx] = :purple
        p.nlabels[][idx] = "Source: "*repr(idx)*"\nNode degree: "*repr(sum(M[idx,:]))
    elseif p.node_color[][idx] == :green
        p.node_color[][idx] = :green1
        p.nlabels[][idx] = "Load: "*repr(idx-num_sources)*"\nNode degree: "*repr(sum(M[idx,:]))
    elseif p.node_color[][idx] == :green1
        p.node_color[][idx] = :green
    end

    p.node_color[] = p.node_color[]
    p.nlabels[] = p.nlabels[]
    p.nlabels_distance[] = p.nlabels_distance[]
end

function edge_click_action(idx, args...)
    if p.edge_color[][idx] == :cyan
        p.edge_color[][idx] = :blue
    else
        p.edge_color[][idx] = :cyan
    end
    p.edge_color[] = p.edge_color[]
end

function node_drag_action(state, idx, event, axis)
    p[:node_pos][][idx] = event.data
    p[:node_pos][] = p[:node_pos][]
end

mutable struct EdgeDragAction
    init::Union{Nothing, Point2f} # save click position
    src::Union{Nothing, Point2f}  # save src vertex position
    dst::Union{Nothing, Point2f}  # save dst vertex position
    EdgeDragAction() = new(nothing, nothing, nothing)
end
function (action::EdgeDragAction)(state, idx, event, axis)
    edge = collect(edges(cm))[idx]
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

function Draw(M, num_sources, num_loads)

    println(num_sources)
    L = num_sources + num_loads

    cm = SimpleGraph(M) # create graph object from Connecitivity Matrix
    
    #_______________________________________________________________________________
    ## defining the colours of the nodes

    nodecolours = Array{Symbol, 1}(undef, L)
    for i in 1:L
        if i <= num_sources
            nodecolours[i] = :red
        else
            nodecolours[i] = :green
        end
    end

    #_______________________________________________________________________________
    ## creating the graph

    edgelabels = Array{String, 1}(undef, ne(cm))
    edgelabels = fill!(edgelabels, "")
    nodelables = Array{String, 1}(undef, nv(cm))
    nodelables = fill!(nodelables, "")

    CM_net, ax, p = graphplot(cm,
                            elabels = edgelabels,
                            nlabels = nodelables,
                            edge_width = [2.0 for i in 1:ne(cm)],
                            edge_color = [:blue for i in 1:ne(cm)],
                            node_size = [10 for i in 1:nv(cm)],
                            node_color = nodecolours);

    #Circular_Layout() # default is Spring graph

    #_______________________________________________________________________________
    ## Creating the Legend

    elem = MarkerElement(color = :blue, marker = 'π', markersize = 15)
    k = Legend(CM_net[1, 2],
        [elem],
        ["I am Legend"])

    #_______________________________________________________________________________
    ## Creating the buttons

    CM_net[2,1] = layout_buttongrid = GridLayout(tellwidth = false)
    layout_buttonlabels = ["Circular Layout", "Spring Layout"]
    layout_buttons = layout_buttongrid[1, 1:2] = [Button(CM_net, label = l) for l in layout_buttonlabels]

    CM_net[0,1] = add_grid = GridLayout(tellwidth = false)

    add_grid[1,1] = add_buttongrid = GridLayout(tellwidth = false)
    add_buttonlabels = ["Add Load", "Add Source", "Add Edge", "Draw Graph"]
    add_buttons = add_buttongrid[1, 1:4] = [Button(CM_net, label = l) for l in add_buttonlabels]

    #_______________________________________________________________________________
    ## Creating the textboxes

    add_grid[2,1] = add_textgrid = GridLayout(tellwidth = false)
    add_textlabels = ["From Node", "To Node"]
    add_text = add_textgrid[1, 1:2] = [Textbox(CM_net, placeholder = l) for l in add_textlabels]

    #_______________________________________________________________________________
    ## Registering all of the actions that can be performed

    ax.aspect = DataAspect()
    deregister_interaction!(ax, :rectanglezoom)
    nhover = NodeHoverHandler(node_hover_action)
    register_interaction!(ax, :nhover, nhover)
    ehover = EdgeHoverHandler(edge_hover_action)
    register_interaction!(ax, :ehover, ehover)
    nclick = NodeClickHandler(node_click_action)
    register_interaction!(ax, :nclick, nclick)
    eclick = EdgeClickHandler(edge_click_action)
    register_interaction!(ax, :eclick, eclick)
    ndrag = NodeDragHandler(node_drag_action)
    register_interaction!(ax, :ndrag, ndrag)
    #register_interaction!(p, :ndrag, ndrag)
    edrag = EdgeDragHandler(EdgeDragAction())
    register_interaction!(ax, :edrag, edrag)

    #_______________________________________________________________________________
    ## Creating the actions

    on(layout_buttons[1].clicks) do clicks;
        Circular_Layout(p, ax)
    end
    on(layout_buttons[2].clicks) do clicks; 
        p.layout = Spring(Ptype=Float32)
        autolimits!(ax)
    end
    on(add_buttons[1].clicks) do clicks; #Add Load
        M = AddNode(M, L+1, L+1) # add at the end, and add one more
        num_loads = num_loads + 1
    
    end
    on(add_buttons[2].clicks) do clicks; #Add Source
        M = AddNode(M, 1, L+1) # add at the beginning, and add one more
        num_sources = num_sources + 1
    end
    on(add_buttons[3].clicks) do clicks; #Add Source
        M = AddNode(M, 1, L+1) # add at the beginning, and add one more
        num_sources = num_sources + 1
    end
    on(add_buttons[4].clicks) do clicks; #Draw Graph
        CM_net, ax, p = Draw(M, num_sources, num_loads)
        display(CM_net)
    end  

    display(CM_net)
    return CM_net, ax, p
end

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#_______________________________________________________________________________
## Connectivity Matrix Parameters

global num_sources = 3
global num_loads = 3

L = num_sources + num_loads

#_______________________________________________________________________________
## Generating Connectivity Matrix

BA = Array{Float64, 2}(undef, L, L)
BA = fill!(BA, 0)

# create scale-free network
BA = Barabasi_Albert(BA, L)

SW = Array{Float64, 2}(undef, L, L)
SW = fill!(SW, 0)

Z = 2 # coordination number, if Z = L, then we have a fully connected graph.
p = 0.4 # the fraction of random connections

# create small world graph
SW = SmallWorld(SW, L, Z, p)

#_______________________________________________________________________________
## Transitioning to general variable names and graph object

M = SW

#_______________________________________________________________________________
## Draw Graph

#Draw(M, num_sources, num_loads)

L = num_sources + num_loads

cm = SimpleGraph(M) # create graph object from Connecitivity Matrix

#_______________________________________________________________________________
## defining the colours of the nodes

nodecolours = Array{Symbol, 1}(undef, L)
for i in 1:L
    if i <= num_sources
        nodecolours[i] = :red
    else
        nodecolours[i] = :green
    end
end

#_______________________________________________________________________________
## creating the graph

edgelabels = Array{String, 1}(undef, ne(cm))
edgelabels = fill!(edgelabels, "")
nodelables = Array{String, 1}(undef, nv(cm))
nodelables = fill!(nodelables, "")

CM_net, ax, p = graphplot(cm,
                        elabels = edgelabels,
                        nlabels = nodelables,
                        edge_width = [2.0 for i in 1:ne(cm)],
                        edge_color = [:blue for i in 1:ne(cm)],
                        node_size = [10 for i in 1:nv(cm)],
                        node_color = nodecolours);


#Circular_Layout() # default is Spring graph

#_______________________________________________________________________________
## Creating the Legend

elem = MarkerElement(color = :blue, marker = 'π', markersize = 15)
k = Legend(CM_net[1, 2],
    [elem],
    ["I am Legend"])

#_______________________________________________________________________________
## Creating the buttons

CM_net[2,1] = layout_buttongrid = GridLayout(tellwidth = false)
layout_buttonlabels = ["Circular Layout", "Spring Layout"]
layout_buttons = layout_buttongrid[1, 1:2] = [Button(CM_net, label = l) for l in layout_buttonlabels]

CM_net[0,1] = add_grid = GridLayout(tellwidth = false)

add_grid[1,1] = add_buttongrid = GridLayout(tellwidth = false)
add_buttonlabels = ["Add Load", "Add Source", "Add Edge", "Draw Graph"]
add_buttons = add_buttongrid[1, 1:4] = [Button(CM_net, label = l) for l in add_buttonlabels]

#_______________________________________________________________________________
## Creating the textboxes

add_grid[2,1] = add_textgrid = GridLayout(tellwidth = false)
add_textlabels = ["From Node", "To Node"]
add_text = add_textgrid[1, 1:2] = [Textbox(CM_net, placeholder = l) for l in add_textlabels]

#_______________________________________________________________________________
## Registering all of the actions that can be performed

ax.aspect = DataAspect()
deregister_interaction!(ax, :rectanglezoom)
nhover = NodeHoverHandler(node_hover_action)
register_interaction!(ax, :nhover, nhover)
ehover = EdgeHoverHandler(edge_hover_action)
register_interaction!(ax, :ehover, ehover)
nclick = NodeClickHandler(node_click_action)
register_interaction!(ax, :nclick, nclick)
eclick = EdgeClickHandler(edge_click_action)
register_interaction!(ax, :eclick, eclick)
ndrag = NodeDragHandler(node_drag_action)
register_interaction!(ax, :ndrag, ndrag)
edrag = EdgeDragHandler(EdgeDragAction())
register_interaction!(ax, :edrag, edrag)

#_______________________________________________________________________________
## Creating the actions

on(layout_buttons[1].clicks) do clicks;
    Circular_Layout(p, ax)
end
on(layout_buttons[2].clicks) do clicks; 
    p.layout = Spring(Ptype=Float32)
    autolimits!(ax)
end
on(add_buttons[1].clicks) do clicks; #Add Load
    M = collect(adjacency_matrix(cm))
    L = size(M,1)
    M = AddNode(M, L+1, L+1) # add at the end, and add one more
    global num_loads = num_loads + 1
end
on(add_buttons[2].clicks) do clicks; #Add Source
    M = collect(adjacency_matrix(cm))
    L = size(M,1)
    M = AddNode(M, 1, L+1) # add at the beginning, and add one more
    global num_sources = num_sources + 1
end
on(add_buttons[3].clicks) do clicks; #Add Edge
    
end
on(add_buttons[4].clicks) do clicks; #Draw Graph
    M = collect(adjacency_matrix(cm))
    Draw(M, num_sources, num_loads)
end  

display(CM_net)

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")
