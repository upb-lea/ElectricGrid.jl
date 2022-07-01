using Stipple, StippleUI, StipplePlotly, PlotlyBase, DataFrames, Graphs, GraphPlot, LinearAlgebra
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
      push!(color_map, "#72C8A9")
    elseif source["fltr"] == "LC"
      push!(color_map, "#E21912")
    elseif source["fltr"] == "L"
      push!(color_map, "#F04842")
    end
  end

  for load in parameters["load"]
    push!(node_descriptions, "Load: " * load["impedance"])

    if load["impedance"] == "RLC"
      push!(color_map, "#400575")
    elseif load["impedance"] == "LC"
      push!(color_map, "#2C128A")
    elseif load["impedance"] == "RL"
      push!(color_map, "#213A95")
    elseif load["impedance"] == "RC"
      push!(color_map, "#35669A")
    elseif load["impedance"] == "L"
      push!(color_map, "#2E99B1")
    elseif load["impedance"] == "C"
      push!(color_map, "#14B6AE")
    elseif load["impedance"] == "R"
      push!(color_map, "#00C979")
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
            width=2
            )
      )
  )

  # Create Plot
  pl = PlotlyBase.Plot(
      [edges_trace, nodes_trace],
      PlotlyBase.Layout(
          hovermode="closest",
          showlegend=false,
          showarrow=false,
          dragmode="select",
          xaxis=attr(showgrid=false, zeroline=false, showticklabels=false),
          yaxis=attr(showgrid=false, zeroline=false, showticklabels=false)
      )
  )

  global model.plot[] = pl
  model.selectedNodes[] = []

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
    push!(parameters["source"], Dict{Any, Any}("L1" => 0.0023, "R_C" => 0.4, "C" => 1.0e-5, "R1" => 0.4, "fltr" => "LC"))
    loc = length(parameters["source"])
  else
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

function deleteNode(pointNumber::Int)
  global CM
  global parameters

  #increment pointNumber since Plotly starts to count at zero
  loc = pointNumber + 1

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

#pl = PlotlyBase.Plot(scatter(x = df.a, text = df.b))

@reactive! mutable struct Setup <: ReactiveModel
  @mixin plot::PBPlotWithEvents()
  showplot::R{Bool} = true
  newSourceTrigger::R{Bool} = false
  newLoadTrigger::R{Bool} = false
  addEdgesTrigger::R{Bool} = false
  deleteNodeTrigger::R{Int64} = 0
  selectedNodes::R{Vector{Dict{String, Any}}} = []
end

#Genie.config.run_as_server = true
Genie.Router.delete!(:Setup)
Stipple.js_mounted(::Setup) = watchplots(Setup)

function handlers(model)

    on(model.isready) do isready
      isready || return
      push!(model)
    end

    on(model.plot_selected) do data
      if haskey(data, "points")
        model.selectedNodes[] = data["points"]
      end
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
            " {{sn.pointIndex}} "
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
  html(ui(model), context = @__MODULE__)
end

up()