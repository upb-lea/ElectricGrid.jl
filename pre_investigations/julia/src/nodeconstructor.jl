using Distributions
using LinearAlgebra
using StatsBase
mutable struct NodeConstructor
    num_connections
    num_source
    num_load
    num_LCL
    num_LC
    num_L
    num_fltr
    cntr
    tot_ele
    CM
    parameter
    S2S_p
    S2L_p
end

function NodeConstructor(;num_source, num_load, CM=nothing, parameter=nothing, S2S_p=0.1, S2L_p=0.8)

    tot_ele = num_source + num_load

    cntr = 0
    num_connections = 0

    if CM === nothing
        cntr, CM = generate_CM(tot_ele, num_source, S2L_p, S2S_p)
        num_connections = cntr
    else
        if size(CM)[1] != tot_ele
            throw("Expect the number of elements in the node to match the specified structure in the CM, but got $tot_ele and $(size(A)[1])")
        end
        num_connections = Int(maximum(CM))
    end

    if parameter === nothing

        #sample = 0.1 * num_source * random.normal(0,1)
        #num_LC = int(np.ceil(np.clip(sample, 1, num_source-1)))

        num_LC = 1 
        num_LCL = num_source - num_LC
        num_L = 0

        parameter = generate_parameter(num_LC, num_LCL, num_L, num_connections, num_load)

    elseif isa(parameter, Dict)
        @assert length(keys(parameter)) == 3 "Expect parameter to have the three entries 'cable', 'load' and 'source' but got $(keys(parameter))"

        @assert length(parameter["source"]) == num_source "Expect the number of sources to match the number of sources in the parameters, but got $num_source and $(length(parameter["source"]))"

        @assert length(parameter["load"]) == num_load "Expect the number of loads to match the number of sources in the parameters, but got $num_load and $(length(parameter["load"]))"

        @assert length(parameter["cable"]) == num_connections "Expect the number of sources to match the number of sources in the parameters, but got $num_connections and $(length(parameter["cable"]))"

        num_LCL, num_LC, num_L = cntr_fltr(parameter["source"])

        @assert num_LCL + num_LC + num_L == num_source "Expect the number of sources to be identical to the sum of the filter types, but the number of sources is $num_source and the sum of the filters is $(num_LCL + num_LC + num_L)"
    else
        throw("Expect parameter to be a dict or nothing, not $(typeof(parameter))")
    end

    num_fltr = 4 * num_LCL + 2 * num_LC + 2 * num_L

    NodeConstructor(num_connections, num_source, num_load, num_LCL, num_LC, num_L, num_fltr,
                cntr, tot_ele, CM, parameter, S2S_p, S2L_p)
end


function generate_parameter(num_LC, num_LCL, num_L, num_connections, num_load)
    """Create parameter dict"""

        source_list = []
        cable_list = []
        load_list = []

        for s in 1:num_LCL
            push!(source_list, sample_LCL_para())
        end
        
        for s in 1:num_LC
            push!(source_list, sample_LC_para())
        end
        
        for s in 1:num_L
            push!(source_list, sample_L_para())
        end

        for c in 1:num_connections
            push!(cable_list, sample_cable_para())
        end

        for l in 1:num_load
            push!(load_list, sample_load_para())
        end

        parameter = Dict()
        parameter["source"] = source_list
        parameter["cable"] = cable_list
        parameter["load"] = load_list
    
        parameter
end


function cntr_fltr(source_list)
    cntr_LCL = 0
    cntr_LC = 0
    cntr_L = 0
    
    for (i, source) in enumerate(source_list)
        if source["fltr"] == "LCL"
            cntr_LCL += 1
        elseif source["fltr"] == "LC"
            cntr_LC += 1
        elseif source["fltr"] == "L"
            cntr_L += 1
        end
    end

    return cntr_LCL, cntr_LC, cntr_L
end

function sample_LCL_para()
    """Sample source parameter""" 

    source = Dict()
    source["fltr"] = "LCL"

    #source["R"] = round(rand(Uniform(0.1, 1)), digits=3)
    #source["L1"] = round(rand(Uniform(2, 2.5)), digits=3) * 1e-3
    #source["L2"] = round(rand(Uniform(2, 2.5)), digits=3) * 1e-3
    #source["C"] = round(rand(Uniform(5, 15)), digits=3) * 1e-6

    #TODO
    source["R"] = 0.4
    source["L1"] = 2.3e-3
    source["L2"] = 2.3e-3
    source["C"] = 10e-6

    source
end

function sample_LC_para()
    """Sample source parameter"""      

    source = Dict()
    source["fltr"] = "LC"

    #source["R"] = round(rand(Uniform(0.1, 1)), digits=3)
    #source["L1"] = round(rand(Uniform(2, 2.5)), digits=3) * 1e-3
    #source["C"] = round(rand(Uniform(5, 15)), digits=3) * 1e-6

    #TODO
    source["R"] = 0.4
    source["L1"] = 2.3e-3
    source["C"] = 10e-6

    source
end

function sample_L_para()
    """Sample source parameter"""      

    source = Dict()
    source["fltr"] = "L"
    source["R"] = round(rand(Uniform(0.1, 1)), digits=3)
    source["L1"] = round(rand(Uniform(2, 2.5)), digits=3) * 1e-3

    source
end


function sample_load_para()
    """Sample load parameter"""

    load = Dict()
    #load["R"] = round(rand(Uniform(10, 10000)), digits=3)

    #TODO
    load["R"] = 14

    load
end

function sample_cable_para()
    """Sample cable parameter"""

    l = rand(1:1:100)

    Rb = 0.722
    Cb = 8*10e-9 # too small?
    Lb = 0.955*10e-3

    cable = Dict()
    cable["R"] = l * Rb
    cable["L"] = l * Lb
    cable["C"] = l * Cb

    #TODO
    cable["R"] = 0.4
    cable["L"] = 2.3e-3
    cable["C"] = 1e-20

    cable
end


function tobe_or_n2b(cntr, x, p)
    """Sets x based on p to zero or to the value of the counter and increments it."""

    # To count up the connection, cntr is returned.
    # If only one type of cable is used this is not necessary an can be replaced by 1
    
    if x < p
        cntr += 1  
        return cntr, x
    else
        x = 0
        return cntr, x
    end
end

function generate_CM(tot_ele, num_source, S2L_p, S2S_p)
    """Constructs the CM

    Returns the constructed CM and the total number of connections.
    """

    # counting the connections 
    cntr = 0

    # get a upper triangular matrix
    mask = UpperTriangular(ones(tot_ele, tot_ele))
    CM = rand(tot_ele,tot_ele) .* mask # fill matrix with random entries between [0,1]
    CM = CM - Diagonal(CM) # delet diagonal bc no connection with itself

    # go through the matrix
    # -1 bc last entry is 0 anyway
    for i in 1:tot_ele-1
        # start at i, bc we need to check only upper triangle
        for j in i:tot_ele-1
            if j >= num_source-1 # select propability according to column
                cntr, x = tobe_or_n2b(cntr, CM[i, j+1], S2L_p)
                CM[i, j+1] = x
            else
                cntr, x = tobe_or_n2b(cntr, CM[i, j+1], S2S_p)
                CM[i, j+1] = x
            end
        end
    end

    # make sure that no objects disappear or subnets are formed
    for i in 1:tot_ele
        # save rows and columns entries
        Col = CM[1:i-1, i]
        Row = CM[i, i+1:tot_ele]
        
        # get one list in the form of: [column, row]-entries
        entries = vcat(Col, Row)

        non_zero = count(i->(i != 0), entries) # number of non_zero entries
        zero = count(i->(i == 0), entries) # number of zero entries

        val_to_set = min(2, zero) # minimum of connections is 2
        
        if non_zero <= 2 # we need to set values if there are less then 2 entries
            idx_row_entries = findall(x->x==0, Col) # Get rows of the entries = 0
            idx_col_entries = findall(x->x==0, Row) # Get col of the entries = 0

            idx_list = vcat([(j,i) for j in idx_row_entries], [(i,i+j) for j in idx_col_entries])
            
            samples = min(val_to_set, length(idx_list))
            idx_rnd = sample(1:length(idx_list), samples) # draw samples from the list
            
            for (a, ix) in enumerate(idx_rnd)
                # Based on the random sample, select an indize
                # from the list and write into the corresponding CM cell.
                cntr += 1
                CM[idx_list[ix][1], idx_list[ix][2]] = cntr
            end
        end
    end

    
    CM = CM - CM' # copy with negative sign to lower triangle

    return cntr, CM
end

function get_A_source(self::NodeConstructor, source_i)
    """Create the A_source entry for a source in the A matrix

    Returns:
        A_source: Matrix with values belonging to corresponding source (2, 2)
    """
    parameter_i = self.parameter["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_source = zeros(4,4)
        A_source[1,1] = -parameter_i["R"]/parameter_i["L1"]
        A_source[1,2] = -1/parameter_i["L1"]
        A_source[2,1] = 1/parameter_i["C"]
        A_source[2,3] = -1/parameter_i["C"]
        A_source[3,2] = 1/parameter_i["L2"]
        A_source[3,4] = -1/parameter_i["L2"]
        
        C_sum =  0
        
        CM_row = self.CM[source_i,:]
        
        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameter["cable"][idx]["C"]
        end
        
        A_source[4,3] = C_sum^(-1)

    elseif parameter_i["fltr"] == "LC"

        A_source = zeros(2,2)
        A_source[1,1] = -parameter_i["R"]/parameter_i["L1"]
        A_source[1,2] = -1/parameter_i["L1"]
        
        C_sum =  parameter_i["C"]
        
        CM_row = self.CM[source_i, :]
        
        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for idx in indizes_
            idx = Int(idx)
            C_sum += self.parameter["cable"][idx]["C"]
        end
        
        A_source[2,1] = C_sum^(-1)

    elseif parameter_i["fltr"] == "L"
        throw("unimplemented")
    else
        throw("Expect filter to be \"LCL\", \"LC\" or \"L\", not $(parameter_i["fltr"]).")
    end

    return A_source

end

function get_B_source(self::NodeConstructor, source_i)
    """Create the B_source entry for a source in the B matrix

    Return:
        B_source: Matrix with values belonging to corresponding source (2, 1)
    """
    parameter_i = self.parameter["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        B_source = zeros(4,1)
        B_source[1,1] =  1/parameter_i["L1"]

    elseif parameter_i["fltr"] == "LC"
        B_source = zeros(2,1)
        B_source[1,1] =  1/parameter_i["L1"]

    elseif parameter_i["fltr"] == "L"
        throw("unimplemented")
    end

    return B_source
end

function get_A_col(self::NodeConstructor, source_i)
    """Create the A_col entry in the A matrix
    Return:
        A_col: Matrix with the column entries for A (2, num_connections)
    """

    parameter_i = self.parameter["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_col = zeros(4, self.num_connections)

        CM_row = self.CM[source_i,:]

        indizes = CM_row[CM_row .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        C_sum = 0 

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameter["cable"][idx]["C"]
        end

        for (idx, sign) in zip(indizes_, signs)
            A_col[4,idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["fltr"] == "LC"

        A_col = zeros(2, self.num_connections)

        CM_row = self.CM[source_i,:]

        indizes = CM_row[CM_row .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices

        C_sum = parameter_i["C"]

        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            C_sum += self.parameter["cable"][idx]["C"]
        end

        for (idx, sign) in zip(indizes_, signs)
            A_col[2,idx] = sign * -(C_sum^(-1))
        end

    elseif parameter_i["fltr"] == "L"
        throw("unimplemented")
    end

    return A_col
end

function get_A_row(self::NodeConstructor, source_i)
    """Create the A_row entry in the A matrix
    Return:
        A_row: Matrix with the row entries for A (num_connections, 2)
    """
    parameter_i = self.parameter["source"][source_i]

    if parameter_i["fltr"] == "LCL"
        A_row = zeros(4, self.num_connections)
        CM_col = self.CM[source_i,:]
        
        indizes = CM_col[CM_col .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_row[4,idx] = sign * 1/self.parameter["cable"][idx]["L"]
        end

    elseif parameter_i["fltr"] == "LC"
        A_row = zeros(2, self.num_connections)
        
        CM_col = self.CM[source_i,:]
        
        indizes = CM_col[CM_col .!= 0] # get entries unequal 0
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        
        for (idx, sign) in zip(indizes_, signs)
            idx = Int(idx)
            A_row[2,idx] = sign * 1/self.parameter["cable"][idx]["L"]
        end

    elseif parameter_i["fltr"] == "L"
        throw("unimplemented")
    end

    return A_row'
end

function generate_A_tran_diag(self::NodeConstructor)
    """Create A_tran_diag"""

    vec = zeros(self.num_connections)
    for (i, ele) in enumerate(self.parameter["cable"])
        vec[i] = -ele["R"]/ele["L"]
    end
    A_tran_diag = Diagonal(vec)

    return A_tran_diag
end

function generate_A_load_col(self::NodeConstructor, load_i)
        
    A_load_col = zeros(self.num_connections)

    CM_row = self.CM[self.num_source + load_i, :]

    indizes = CM_row[CM_row .!= 0] # get entries unequal 0
    signs = [sign(x) for x in indizes] # get signs

    indizes_ = indizes .* signs # delet signs from indices

    C_sum = 0

    for ele in self.parameter["cable"]
        C_sum += ele["C"]
    end

    for (idx, sign) in zip(indizes_, signs)
        idx = Int(idx)
        A_load_col[idx] = sign * -(C_sum^(-1))
    end

    return A_load_col
end


function generate_A_load_row(self::NodeConstructor, load_i)
            
    A_load_row = zeros(self.num_connections)

    CM_col = self.CM[self.num_source + load_i, :]

    indizes = CM_col[CM_col .!= 0] # get entries unequal 0
    signs = [sign(x) for x in indizes] # get signs
    indizes_ = indizes .* signs # delet signs from indices

    for (idx, sign) in zip(indizes_, signs)
        idx = Int(idx)

        A_load_row[idx] = sign * 1/self.parameter["cable"][idx]["L"] 
    end 

    return A_load_row
end


function generate_A_load_diag(self::NodeConstructor)
    vec = zeros(self.num_load)

    for (i, ele) in enumerate(self.parameter["load"])
        CM_row = self.CM[self.num_source + i, :]
        indizes = CM_row[CM_row .!= 0]
        signs = [sign(x) for x in indizes] # get signs
        indizes_ = indizes .* signs # delet signs from indices
        C_sum = 0

        for idx in indizes_
            idx = Int(idx)
            
            C_sum += self.parameter["cable"][idx]["C"]  # Cb[idx]
        end

        vec[i] = - (ele["R"] * (C_sum))^(-1) # Rload[i]
    end

    A_load_diag = Diagonal(vec)

    return A_load_diag
end


function generate_A(self::NodeConstructor)
    """Generate the A matrix

    The previously constructed matrices are now plugged together in the form:
        [[A_source, A_col,       0          ],
        [A_row,    A_tran_diag, A_load_row ],
        [0,        A_load_col,  A_load_diag]]
    with A_source:
        [[LCL, 0 , 0]
        [0,   LC, 0]
        [0,   0,  L]]
    Returns:
        A: A matrix for state space ((2*num_source+num_connections),(2*num_source+num_connections))
    """
    # get A_source
    self.num_fltr = 4 * self.num_LCL + 2 * self.num_LC + 2 * self.num_L
    A_source = zeros(self.num_fltr, self.num_fltr) # construct matrix of zeros
    A_source_list = [get_A_source(self, i) for i in 1:self.num_source]
            
    for (i, ele) in enumerate(A_source_list)
        if i <= self.num_LCL
            start = 4 * (i - 1) + 1
            stop = 4 * i
            A_source[start:stop, start:stop] = ele

        elseif i <= self.num_LCL + self.num_LC
            start = 2 * i + 2 * self.num_LCL - 1
            stop = 2 * i + 2 * self.num_LCL
            A_source[start:stop, start:stop] = ele

        elseif i <= self.num_LCL + self.num_LC + self.num_L
            start = 2 * i + 2 * self.num_LCL - 1
            stop = 2 * i + 2 * self.num_LCL
            A_source[start:stop,start:stop] = ele
        end
    end

    # get A_col
    A_col = zeros(self.num_fltr, self.num_connections)
    A_col_list = [get_A_col(self, i) for i in 1:self.num_source] # start at 1 bc Source 1 ...

    for (i, ele) in enumerate(A_col_list)
        if i <= self.num_LCL
            start = 4 * (i - 1) + 1
            stop = 4 * i
            A_col[start:stop,:] = ele

        elseif i <= self.num_LCL+self.num_LC
            start = 2 * i + 2 * self.num_LCL - 1
            stop = 2 * i + 2 * self.num_LCL
            A_col[start:stop,:] = ele

        elseif i <= self.num_LCL+self.num_LC+self.num_L
            start = 2 * i + 2 * self.num_LCL - 1
            stop = 2 * i + 2 * self.num_LCL
            A_col[start:stop,:] = ele
        end
    end

    # get A_row
    A_row = zeros(self.num_connections, self.num_fltr)
    A_row_list = [get_A_row(self, i) for i in 1:self.num_source] # start at 1 bc Source 1 ...

    for (i, ele) in enumerate(A_row_list)
        if i <= self.num_LCL
            start = 4 * (i - 1) + 1
            stop = 4 * i
            A_row[:,start:stop] = ele

        elseif i <= self.num_LCL+self.num_LC
            start = 2 * i + 2 * self.num_LCL - 1
            stop = 2 * i + 2 * self.num_LCL
            A_row[:,start:stop] = ele

        elseif i <= self.num_LCL+self.num_LC+self.num_L
            start = 2 * i + 2 * self.num_LCL - 1
            stop = 2 * i + 2 * self.num_LCL
            A_row[:,start:stop] = ele
        end
    end

    A_tran_diag = generate_A_tran_diag(self)

    A_load_row_list = []
    for i in 1:self.num_load
        push!(A_load_row_list, generate_A_load_row(self, i))
    end
    A_load_row = reduce(hcat, A_load_row_list) # i-> idx // i+1 -> num of load

    A_load_col_list = []
    for i in 1:self.num_load
        push!(A_load_col_list, generate_A_load_col(self, i))
    end
    A_load_col = reduce(hcat, A_load_col_list)'

    A_load_diag = generate_A_load_diag(self)

    #         A_transition = np.block([A_tran_diag, A_load_col],
    #                                 [A_load_col, A_load_diag])

    A_load_zeros = zeros(self.num_fltr, self.num_load)

    A_load_zeros_t = A_load_zeros'

    A = [A_source       A_col        A_load_zeros
        A_row           A_tran_diag  A_load_row
        A_load_zeros_t  A_load_col   A_load_diag]

    return A
end



function generate_B(self::NodeConstructor)
    """Generate the B matrix

    The previously constructed matrices are now plugged together in the form:
        [[B_source,        0, ...,         0],
        [       0, B_source, ...,         0],
        [       0,        0, ...,         0],
        [       0,        0, ...,  B_source]]
    Returns:
        B: B matrix for state space (2*num_source+num_connections,num_source)
    """
    B = zeros(self.num_fltr + self.num_connections + self.num_load, self.num_source)

    B_source_list = [get_B_source(self, i) for i in 1:self.num_source] # start at 1 bc Source 1 ...

    for (i, ele) in enumerate(B_source_list)
        if i <= self.num_LCL
            start = start = 4 * (i - 1) + 1
            stop = 4 * i
            B[start:stop,i] = ele

        elseif i <= self.num_LCL + self.num_LC
            start = 2 * i + 2 * self.num_LCL - 1
            stop = 2 * i + 2 * self.num_LCL
            B[start:stop,i] = ele

        elseif i <= self.num_LCL + self.num_LC + self.num_L
            start = 2 * i + 2 * self.num_LCL - 1
            stop = 2 * i + 2 * self.num_LCL
            B[start:stop,i:i] = ele
        end
    end

    return B
end



function generate_C(self::NodeConstructor)
    """Generate the C matrix

    Retruns:
        C: Identity matrix (2*num_source+num_connections)
    """
    return Diagonal(ones(self.num_fltr + self.num_connections + self.num_load))
end

function generate_D(self::NodeConstructor)
    """Generate the D vector

    Retruns:
        0: Zero vector (2*num_source+num_connections)
    """
    return 0
end

function get_sys(self::NodeConstructor)
    """Returns state space matrices"""

    A = generate_A(self)
    B = generate_B(self)
    C = generate_C(self)
    D = generate_D(self)
    return (A, B, C, D)
end


function get_states(self::NodeConstructor)
    states = []
    for s in 1:self.num_source
        if s <= self.num_LCL
            push!(states, "i_f$s")    # i_f1; dann i_f2....
            push!(states, "u_f$s")
            push!(states, "i_$s")
            push!(states, "u_$s")

        elseif s <= self.num_LCL + self.num_LC + self.num_L
            push!(states, "i_$s")
            push!(states, "u_$s")
        end
    end

    for c in 1:self.num_connections
        push!(states, "i_c$c")
    end

    for l in 1:self.num_load
        push!(states, "u_l$l")
    end
    return states
end

function draw_graph(self::NodeConstructor)
    """Plots a graph according to the CM matrix

    Red nodes corresponse to a source.
    Lightblue nodes corresponse to a load.
    """

    # edges = []
    # color = []
    # for i in range(1, self.num_connections+1):
    #     (row, col) = np.where(self.CM==i)
    #     (row_idx, col_idx) = (row[0]+1, col[0]+1)
    #     edges.append((row_idx, col_idx))
    #     if row_idx <= self.num_source:
    #         color.append('red')
    #     else:
    #         color.append('blue')
    #     end
    # end

    # G = nx.Graph(edges)

    # color_map = []

    # for node in G:
    #     if node <= self.num_source:
    #         color_map.append("red")
    #     else:
    #         color_map.append("lightblue")
    #     end
    # end

    # nx.draw(G, node_color=color_map, with_labels = True)
    # plt.show()

    # pass
end
