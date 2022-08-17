using Stipple, StippleUI, StipplePlotly, PlotlyBase, DataFrames, Graphs, GraphPlot, LinearAlgebra, GenieAutoReload
import Genie.Renderer.Html.div

register_mixin(@__MODULE__)

global model
global CM
global parameters

function drawGraph()
  global CM
  global parameters

  CMtemp = CM + -2 * LowerTriangular(CM)

  G = SimpleGraph(CMtemp)

  # Position nodes
  pos_x, pos_y = GraphPlot.shell_layout(G)

  # Create plot points
  edge_x = []
  edge_y = []

  for edge in edges(G)
      push!(edge_x, pos_x[src(edge)])
      push!(edge_x, pos_x[dst(edge)])
      push!(edge_x, nothing)
      push!(edge_y, pos_y[src(edge)])
      push!(edge_y, pos_y[dst(edge)])
      push!(edge_y, nothing)
  end

  #  Color nodes
  color_map = []
  node_descriptions = []

  for source in parameters["source"]
    push!(node_descriptions, "Source: " * source["fltr"])

    if source["fltr"] == "LCL"
      push!(color_map, "#FF8800")
    elseif source["fltr"] == "LC"
      push!(color_map, "#FF6600")
    elseif source["fltr"] == "L"
      push!(color_map, "#FF3300")
    end
  end

  for load in parameters["load"]
    push!(node_descriptions, "Load: " * load["impedance"])

    if load["impedance"] == "RLC"
      push!(color_map, "#8F00D1")
    elseif load["impedance"] == "LC"
      push!(color_map, "#4900A8")
    elseif load["impedance"] == "RL"
      push!(color_map, "#3A09C0")
    elseif load["impedance"] == "RC"
      push!(color_map, "#0026FF")
    elseif load["impedance"] == "L"
      push!(color_map, "#0066FF")
    elseif load["impedance"] == "C"
      push!(color_map, "#00CCFF")
    elseif load["impedance"] == "R"
      push!(color_map, "#00F3E7")
    end
  end
  

  # Create edges
  edges_trace = scatter(
      mode="lines",
      x=edge_x,
      y=edge_y,
      line=attr(
          width=0.8,
          color="#113"
      ),
  )

  # Create nodes
  nodes_trace = scatter(
      x=pos_x,
      y=pos_y,
      mode="markers",
      text = node_descriptions,
      marker=attr(
          color=color_map,
          size=13,
          line=attr(
            color="Black",
            width=1
            )
      )
  )

  # Create Plot
  pl = PlotlyBase.Plot(
      [edges_trace, nodes_trace],
      PlotlyBase.Layout(
          plot_bgcolor="#f1f3f7",
          hovermode="closest",
          showlegend=false,
          showarrow=false,
          dragmode="select",
          xaxis=attr(showgrid=false, zeroline=false, showticklabels=false),
          yaxis=attr(showgrid=false, zeroline=false, showticklabels=false)
      )
  )

  global model.plot[] = pl
  #model.selectedNodes[] = []

  nothing
end

function startConfig()
  global CM = [ 0   1
       -1   0];

  # Fill parameter dict
  cable_list = []
  load_list = []
  source_list = []

  # LC filter
  source = Dict()
  source["fltr"] = "LC"
  source["R1"] = 0.4
  #source["R2"] = 0.4
  source["R_C"] = 0.4
  source["L1"] = 2.3e-3
  #source["L2"] = 2.3e-3
  source["C"] = 10e-6;

  # Load
  load = Dict()
  load["impedance"] = "R"
  load["R"] = 14;

  # cable 
  cable = Dict()
  cable["R"] = 0.722
  cable["L"] = 0.264e-3
  cable["C"] = 0.4e-6


  push!(source_list, source);
  push!(load_list, load);
  push!(cable_list, cable);

  global parameters = Dict()
  parameters["source"] = source_list
  parameters["cable"] = cable_list
  parameters["load"] = load_list

  nothing
end

function addNode(type::String)
  global CM
  global parameters

  if type == "source"
    # loc = 0
    # for i = 1:length(parameters["source"])
    #   if parameters["source"][i]["fltr"] != "L"
    #     loc = i
    #     break
    #   end
    # end
    # if loc == 0
    #   loc = length(parameters["source"])
    # end
    # insert!(parameters["source"], loc, Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC"))
    push!(parameters["source"], Dict{Any, Any}("pwr" => 10000.0, "vdc" => 700, "i_rip" => 0.12, "v_rip" => 0.015, "L1" => 0.0023, "R1" => 0.4, "fltr" => "L"))
    loc = length(parameters["source"])
  else
    # loc = 0
    # for i = 1:length(parameters["load"])
    #   if parameters["source"][i]["impedance"] != "L"
    #     loc = length(parameters["source"]) + i
    #     break
    #   end
    # end
    # if loc == 0
    #   loc = length(parameters["source"]) + length(parameters["load"])
    # end
    push!(parameters["load"], Dict{Any, Any}("R" => 14, "impedance" => "R"))
    loc = length(parameters["source"]) + length(parameters["load"])
  end

  L_or = size(CM,1)
  L_new = L_or + 1

  new_M = zeros(Int64, L_new, L_new)

  new_M[1:loc-1,1:loc-1] = CM[1:loc-1,1:loc-1]
  new_M[loc+1:end,1:loc-1] = CM[loc:end,1:loc-1]
  new_M[1:loc-1,loc+1:end] = CM[1:loc-1,loc:end]
  new_M[loc+1:end,loc+1:end] = CM[loc:end,loc:end]

  global CM = new_M

  drawGraph()
  nothing
end

function addEdges()
  #always add +1 to pointIndex since Plotly starts to count at zero
  global CM
  global parameters
  global model

  nextnumber = maximum(CM) + 1

  #flip CM to all-positive
  CMtemp = CM + -2 * LowerTriangular(CM)

  for pair in Iterators.product(model.selectedNodes[], model.selectedNodes[])
    #check if both nodes aren't the same
    if pair[1]["pointIndex"] != pair[2]["pointIndex"]
      #check if no edge exists
      if CMtemp[pair[1]["pointIndex"]+1, pair[2]["pointIndex"]+1] == 0
        CMtemp[pair[1]["pointIndex"]+1, pair[2]["pointIndex"]+1] = nextnumber
        CMtemp[pair[2]["pointIndex"]+1, pair[1]["pointIndex"]+1] = nextnumber
        push!(parameters["cable"], Dict{Any, Any}("C" => 4.0e-7, "L" => 0.000264, "R" => 0.722))
        nextnumber += 1
      end
    end
  end

  #flip back
  global CM = CMtemp + -2 * LowerTriangular(CMtemp)

  drawGraph()
  nothing
end

function deleteNode(pointIndex::Int)
  global CM
  global parameters

  loc = pointIndex + 1
  i = 1
  selectedNodes = model.selectedNodes[]
  while i <= length(selectedNodes)
    found = false
    if selectedNodes[i]["pointIndex"] == pointIndex
      deleteat!(selectedNodes, i)
      found = true
    elseif found
      selectedNodes[i]["pointIndex"] -= 1
      selectedNodes[i]["pointNumber"] -= 1
    end
    i += 1
  end
  model.selectedNodes[] = selectedNodes
  #model.selectedNodes[] = []

  if loc <= length(parameters["source"])
    deleteat!(parameters["source"], loc)
  else
    deleteat!(parameters["load"], loc - length(parameters["source"]))
  end

  deletedEdges = filter( (x) -> x != 0, CM[:,loc] )

  L_or = size(CM,1)
  L_new = L_or - 1
  new_M = zeros(Int64, L_new, L_new)

  new_M[1:loc-1,1:loc-1] = CM[1:loc-1,1:loc-1]
  new_M[loc:end,1:loc-1] = CM[loc+1:end,1:loc-1]
  new_M[1:loc-1,loc:end] = CM[1:loc-1,loc+1:end]
  new_M[loc:end,loc:end] = CM[loc+1:end,loc+1:end]

  deletedEdges = abs.(deletedEdges)
  sort!(deletedEdges)
  deletedEdgesReverse = sort(deletedEdges, rev=true)

  for ele in deletedEdgesReverse
    deleteat!(parameters["cable"], ele)
  end

  #fix edge numbers
  for (i, ele) in enumerate(deletedEdges)
    f(x) = x > (ele + 1 - i) ? x - 1 : x
    g(x) = x < - (ele + 1 - i) ? x + 1 : x
    new_M = f.(new_M)
    new_M = g.(new_M)
  end

  global CM = new_M

  drawGraph()
  nothing
end

function swapcols!(X::AbstractMatrix, i::Integer, j::Integer)
  @inbounds for k = 1:size(X,1)
      X[k,i], X[k,j] = X[k,j], X[k,i]
  end
end

function swaprows!(X::AbstractMatrix, i::Integer, j::Integer)
  @inbounds for k = 1:size(X,2)
      X[i,k], X[j,k] = X[j,k], X[i,k]
  end
end

#pl = PlotlyBase.Plot(scatter(x = df.a, text = df.b))

@reactive! mutable struct Setup <: ReactiveModel
  @mixin plot::PBPlotWithEvents()
  showplot::R{Bool} = true
  newSourceTrigger::R{Bool} = false
  newLoadTrigger::R{Bool} = false
  addEdgesTrigger::R{Bool} = false
  deleteNodeTrigger::R{Int64} = 0
  selectedNodes::R{Vector{Dict{String, Any}}} = []
  sourcemodels::R{Vector{String}} = ["LCL", "LC", "L"]
  loadmodels::R{Vector{String}} = ["RLC", "LC", "RL", "L", "RC", "C", "R"]
  test::R{String} = "LC"
end

#Genie.config.run_as_server = true
Genie.Router.delete!(:Setup)
Stipple.js_mounted(::Setup) = watchplots()

function handlers(model)

    on(model.isready) do isready
      isready || return
      push!(model)
    end

    on(model.plot_selected) do data
      if haskey(data, "points")
        selectedNodes = data["points"]

        for i = 1:length(selectedNodes)
          tempnode = selectedNodes[i]
          if tempnode["pointIndex"] + 1 <= length(parameters["source"])
            tempnode["type"] = "source"
            tempnode["parameters"] = parameters["source"][tempnode["pointIndex"] + 1]
          else
            tempnode["type"] = "load"
            tempnode["parameters"] = parameters["load"][tempnode["pointIndex"] + 1 - length(parameters["source"])]
          end
          selectedNodes[i] = tempnode
        end

        model.selectedNodes[] = selectedNodes
      end
    end

    on(model.selectedNodes) do _
      println("node change detected")

      #TODO col and row switch, also what happens to the selectedNodes?
      #-> Alle SelectedNodes von hinten durchgehen und aus den Dicts rauslöschen, dann neu einordnen

      for i = 1:length(model.selectedNodes[])
        tempnode = model.selectedNodes[i]
        if tempnode["pointIndex"] + 1 <= length(parameters["source"])
          parameters["source"][tempnode["pointIndex"] + 1] = model.selectedNodes[i]["parameters"]
        else
          parameters["load"][tempnode["pointIndex"] + 1 - length(parameters["source"])] = model.selectedNodes[i]["parameters"]
        end
      end

      drawGraph()
    end

    on(model.newSourceTrigger) do _
      if (model.newSourceTrigger[])
        addNode("source")
        model.newSourceTrigger[] = false
      end
    end

    on(model.newLoadTrigger) do _
      if (model.newLoadTrigger[])
        addNode("load")
        model.newLoadTrigger[] = false
      end
    end

    on(model.addEdgesTrigger) do _
      if (model.addEdgesTrigger[])
        addEdges()
        model.addEdgesTrigger[] = false
      end
    end

    on(model.deleteNodeTrigger) do _
      println("DELETED NODE TRIGGER DETECTED!!!!!!!!!!!!!!")
      deleteNode(model.deleteNodeTrigger[])
    end
  
    model
  end

function ui(model)
  page( model, class="container q-layout", [
      
      header(class="st-header q-pa-sm", row([
        Stipple.image(src="logo.png")

        toggle("Microgrid", :showplot)
      ])) 


      row([

        cell(class="st-module", size=9, style="height:600px;", [
          h4("Microgrid Graph")

          btn("Add Source", @click("newSourceTrigger = true"))
          btn("Add Load", @click("newLoadTrigger = true"))
          btn("Add Edges", @click("addEdgesTrigger = true"))

          plotly(:plot, syncevents = true)
        ])

        scrollarea(class="col-2", style="height:600px; overflow: auto;", [
          cell([
            Stipple.p("{{sn.type}}: {{sn.parameters.fltr}}", @iif("sn.type == \"source\""))
            Stipple.p("{{sn.type}}: {{sn.parameters.impedance}}", @iif("sn.type == \"load\""))
            Stipple.p(StippleUI.Selects.select(:var"sn.parameters.fltr", options=:sourcemodels , label="Model"), @iif("sn.type == \"source\""))
            Stipple.p(StippleUI.Selects.select(:var"sn.parameters.impedance", options=:loadmodels , label="Model"), @iif("sn.type == \"load\""))

            Stipple.p([
              textfield("L1", :var"sn.parameters.L1", @iif("sn.parameters.fltr == \"LCL\" || sn.parameters.fltr == \"LC\" || sn.parameters.fltr == \"L\""), rules = "[val => !(isNaN(val))]")
              textfield("L2", :var"sn.parameters.L2", @iif("sn.parameters.fltr == \"LCL\""), rules = "[val => !(isNaN(val))]")
              textfield("C", :var"sn.parameters.C", @iif("sn.parameters.fltr == \"LCL\" || sn.parameters.fltr == \"LC\""), rules = "[val => !(isNaN(val))]")
              textfield("R1", :var"sn.parameters.R1", @iif("sn.parameters.fltr == \"LCL\" || sn.parameters.fltr == \"LC\" || sn.parameters.fltr == \"L\""), rules = "[val => !(isNaN(val))]")
              textfield("R2", :var"sn.parameters.R2", @iif("sn.parameters.fltr == \"LCL\""), rules = "[val => !(isNaN(val))]")
              textfield("R_C", :var"sn.parameters.R_C", @iif("sn.parameters.fltr == \"LCL\" || sn.parameters.fltr == \"LC\""), rules = "[val => !(isNaN(val))]")
            ], @iif("sn.type == \"source\""))

            Stipple.p([
              textfield("R", :var"sn.parameters.R", @iif("sn.parameters.impedance == \"RLC\" || sn.parameters.impedance == \"RL\" || sn.parameters.impedance == \"RC\" || sn.parameters.impedance == \"R\""), rules = "[val => !(isNaN(val))]")
              textfield("L", :var"sn.parameters.L", @iif("sn.parameters.impedance == \"RLC\" || sn.parameters.impedance == \"LC\" || sn.parameters.impedance == \"RL\" || sn.parameters.impedance == \"L\""), rules = "[val => !(isNaN(val))]")
              textfield("C", :var"sn.parameters.C", @iif("sn.parameters.impedance == \"RLC\" || sn.parameters.impedance == \"LC\" || sn.parameters.impedance == \"RC\" || sn.parameters.impedance == \"C\""), rules = "[val => !(isNaN(val))]")
            ], @iif("sn.type == \"load\""))

            btn("delete", @click("deleteNodeTrigger = sn.pointIndex"))
          ], class="st-module", style="margin-bottom:20px; padding:8px;", @recur(:"sn in selectedNodes"))
        ])

      ], class="flex-center", @iif("showplot"))

    ]
  )
end



route("/") do
  global model = Setup |> init |> handlers
  global CM
  global parameters
  if !(@isdefined CM) || !(@isdefined parameters)
    startConfig()
  end
  drawGraph()
  #Stipple.injectdeps(GenieAutoReload.assets(), model)
  html(ui(model), context = @__MODULE__)
end

#Genie.config.websockets_server = true
GenieAutoReload.autoreload(pwd())
up()