using GraphMakie
using Graphs
using GLMakie
using Images

function AddNode(M, node)

    # Adds new unconnected nodes to the system (if it is not already there)

    if size(M,1) < node
        increase = node-size(M,1)
        L = size(M,1)+ increase
        new_M = Array{Float64, 2}(undef, L, L)
        new_M = fill!(new_M, 0)
        new_M[1:L-increase, 1:L-increase] = M
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
    he Barabasi-Albert model is an algorithm for generating random scale-free networks
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
    node. Therefore, we are returnign a vector.
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

function node_hover_action(state, idx, event, axis)
    p.node_size[][idx] = state ? 20 : 10
    p.node_size[] = p.node_size[] # trigger observable
end

function edge_hover_action(state, idx, event, axis)
    p.edge_width[][idx]= state ? 10.0 : 2.0
    p.edge_width[]  = p.edge_width[]  # trigger observable
    
    p.elabels[][idx] = state ? "a" : repr(idx) # ? is the ternary operator, closely related to the if-elseif-else
    p.elabels[] = p.elabels[]
    #p.background[] = p.background[]

    println(p.elabels)
    #display(sw_net)
end

function node_click_action(idx, args...)
    p.node_color[][idx] = rand(RGB)
    p.node_color[] = p.node_color[]
end

function edge_click_action(idx, args...)
    p.edge_color[][idx] = rand(RGB)
    p.edge_color[] = p.edge_color[]
end

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n")

#_______________________________________________________________________________
num_sources = 5
num_loads = 10

L = num_sources + num_loads # The L nodes in a small world networks are arranged around a circle

BA = Array{Float64, 2}(undef, L, L)
BA = fill!(BA, 0)

# create scale-free network
BA = Barabasi_Albert(BA, L)

SW = Array{Float64, 2}(undef, L, L)
SW = fill!(SW, 0)

Z = 2 # coordination number, if Z = L, then we have a fully connected graph.

p = 0.2 # the fraction of random connections

# create small world graph
SW = SmallWorld(SW, L, Z, p)

sw = SimpleGraph(SW)

# define some colors
nodecolours = [:blue for i in 1:nv(sw)]
nodecolours = [:yellow for i in 1:nv(sw)]
nodecolours[4] = nodecolours[7] = :red

for i in 1:L
    if i <= num_sources
        nodecolours[i] = :red
    else
        nodecolours[i] = :green
    end
end

edgelabels = repr.(1:ne(sw))
#edgelabels[1] = "test"
#edgelabels = ["edge $i" for i in 1:ne(sw)]

sw_net, ax, p = graphplot(sw,
                        elabels = edgelabels,
                        edge_width = [2.0 for i in 1:ne(sw)],
                        edge_color = [:blue for i in 1:ne(sw)],
                        node_size = [10 for i in 1:nv(sw)],
                        node_color = nodecolours)

display(sw_net)

#p.elabels = (repr.(200*(1:ne(sw))))
#p.elabels[][1] = "aesdfasdfa "
#println(p.elabels)

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

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")