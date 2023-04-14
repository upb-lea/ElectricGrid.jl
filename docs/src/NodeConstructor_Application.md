# NodeConstructor Application

The creation of a grid always follows a certain pattern and can therefore be automated. In this way, you can avoid errors that may arise from the handwritten rearrangement of the equations. It is therefore obvious to write an automatism for the generation, which generates the matrices on the basis of input parameters. In this notebook we will discuss:

- ### Automatic Generation of Grids
- ### Three-phase Systems

## Automatic Generation of Grids

In this part of the notebook the example of the `NodeConstructor_Theory.ipynb` shall be reproduced with the help of the NodeConstructor.

An important point is the information about how the components are connected by which cable. To pass this information to the NodeConstructor we introduce the connection matrix (`CM`).

### Structured Grid Generation via Connection Matrix

The connection matrix specifies how the elements of the grid are connected. In our grid there are two basic elements: sources and loads, which are then connected via cables.

To get a better understanding of the CM, the CM of the previous example is shown below:

| From\To  | Source 1 | Source 2 | Load 1 |
|:---------|:--------:|:--------:|:------:|
| **Source 1** |     0    |     0    |    1   |
| **Source 2** |     0    |     0    |    2   |
| **Load 1**   |    -1    |    -2    |    0   |

The first column of the above matrix indicates which element we start from. The sources are considered first and after that the loads. The first row of the matrix indicates the elements to which we can connect. The number and order is identical to the first column.

Since our grid elements cannot be connected to themselves, there are always zeros on the main diagonal. In the second row of the table we can now identify the connections of the first source. This source is not connected to the second source but to the first load, which is why we find a "1" in this cell. We can read this entry as: "Source 1 is connected to load 1 by cable 1". In the third row of the table above, all connections are made from source 2, which is connected to load 1 via cable 2. The fourth line contains the connections of the 1 load. As shown before, it is connected to source 1 and source 2, but the entries have a negative sign. The negative sign indicates that the flow direction of the current is assumed to be negative. The sign of the entries is taken into account when constructing the DGLs, but has no influence on the subsequent current flow in the simulation. The matrix is antisymmetric, which means that we have a mirroring of the matrix along the main diagonal with the sign of the entries reversed.

In Julia, this matrix then looks like this:


```julia
CM = [ 0  0  1
       0  0  2
      -1 -2  0];
```

We now import the ElectricGrid package:


```julia
# Import the package
using ElectricGrid
```


```@raw html
<div style="padding: 1em; background-color: #f8d6da; border: 1px solid #f5c6cb; font-weight: bold;">
<p>The WebIO Jupyter extension was not detected. See the
<a href="https://juliagizmos.github.io/WebIO.jl/latest/providers/ijulia/" target="_blank">
    WebIO Jupyter integration documentation
</a>
for more information.
</div>
```



To reproduce the above example exactly, the important parameters must be passed through the parameter dict: 


```julia
# Source
R = 1.1e-3
L = 70e-6
R_c = 7e-3
C = 250e-6

# Cable
C_b = 1e-4/2
L_b = 1e-4
R_b = 1e-3

# Load
R_l = 100
C_l = 1e-2
L_l = 1e-2;
```


```julia
parameters = Dict()

grid_properties = Dict()
grid_properties["fs"] =  10e3
grid_properties["v_rms"] = 230
grid_properties["phase"] = 1;
parameters["grid"] = grid_properties

source1 = Dict()
source2 = Dict()
source_list = []

source1["fltr"] = "LCL"
source1["R1"] = R
source1["L1"] = L
source1["C"] = C
source1["R_C"] = R_c
source1["R2"] = R
source1["L2"] = L

source2["fltr"] = "LC"
source2["R1"] = R
source2["L1"] = L
source2["C"] = C
source2["R_C"] = R_c
push!(source_list, source1, source2)

parameters["source"] = source_list

cable = Dict()
cable["R"] = R_b
cable["L"] = L_b
cable["C"] = C_b
cable_list = []

push!(cable_list, cable, cable);
parameters["cable"] = cable_list

load1 = Dict()
load_list = []

load1["impedance"] = "RLC"
load1["R"] = R_l;
load1["L"] = L_l;
load1["C"] = C_l;

push!(load_list, load1);
parameters["load"] = load_list;
```

Now the NodeConstructor is called. In addition to the number of sources and loads, this also receives the CM and the parameter dict:


```julia
S2_L1 = NodeConstructor(num_sources=2, num_loads=1, parameters=parameters, CM=CM);
```

With the function `DrawGraph()` the topology of the grid can now be displayed. Here, the color orange corresponds to a source and the color blue corresponds to a load.


```julia
DrawGraph(S2_L1)
```


```@raw html
<html>
<head><meta charset="utf-8" /></head>
<body>
<div>
    <script src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/MathJax.js?config=TeX-AMS-MML_SVG"></script>
    <script type="text/javascript">
window.PlotlyConfig = {MathJaxConfig: 'local'};
</script>
<script type="text/javascript">
if (window.MathJax) {MathJax.Hub.Config({SVG: {font: "STIX-Web"}});}
</script>
<script type="text/javascript">
    if (typeof require !== 'undefined') {
        require.undef("plotly");
        requirejs.config({
            paths: {
                'plotly': ['https://cdn.plot.ly/plotly-2.3.0.min']
            }
        });
        require(['plotly'], function(Plotly) {
            window._Plotly = Plotly;
        });
    }
</script>

    <div
        id=98fa629d-3d68-41da-a6b8-33bde237e67a
        class="plotly-graph-div"
        style="height:100%; width:100%;">
    </div>
    <script type="text/javascript">
        require(["plotly"], function(Plotly) {
        window.PLOTLYENV = window.PLOTLYENV || {}

        if (document.getElementById('98fa629d-3d68-41da-a6b8-33bde237e67a')) {
    Plotly.newPlot(
        '98fa629d-3d68-41da-a6b8-33bde237e67a',
        [{"mode":"lines","line":{"color":"#113","width":0.8},"y":[0.0,-0.8660254037844385,null,0.8660254037844387,-0.8660254037844385,null],"type":"scatter","x":[1.0,-0.5000000000000004,null,-0.4999999999999998,-0.5000000000000004,null]},{"x":[1.0,-0.4999999999999998,-0.5000000000000004],"mode":"markers","y":[0.0,0.8660254037844387,-0.8660254037844385],"type":"scatter","text":["Source: LCL","Source: LC","Load: RLC"],"marker":{"color":["#FF8800","#FF6600","#8F00D1"],"line":{"color":"Black","width":1},"size":13}}],
        {"showlegend":false,"xaxis":{"showgrid":false,"showticklabels":false,"zeroline":false},"hovermode":"closest","template":{"layout":{"coloraxis":{"colorbar":{"ticks":"","outlinewidth":0}},"xaxis":{"gridcolor":"white","zerolinewidth":2,"title":{"standoff":15},"ticks":"","zerolinecolor":"white","automargin":true,"linecolor":"white"},"hovermode":"closest","paper_bgcolor":"white","geo":{"showlakes":true,"showland":true,"landcolor":"#E5ECF6","bgcolor":"white","subunitcolor":"white","lakecolor":"white"},"colorscale":{"sequential":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]],"diverging":[[0,"#8e0152"],[0.1,"#c51b7d"],[0.2,"#de77ae"],[0.3,"#f1b6da"],[0.4,"#fde0ef"],[0.5,"#f7f7f7"],[0.6,"#e6f5d0"],[0.7,"#b8e186"],[0.8,"#7fbc41"],[0.9,"#4d9221"],[1,"#276419"]],"sequentialminus":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]},"yaxis":{"gridcolor":"white","zerolinewidth":2,"title":{"standoff":15},"ticks":"","zerolinecolor":"white","automargin":true,"linecolor":"white"},"shapedefaults":{"line":{"color":"#2a3f5f"}},"hoverlabel":{"align":"left"},"mapbox":{"style":"light"},"polar":{"angularaxis":{"gridcolor":"white","ticks":"","linecolor":"white"},"bgcolor":"#E5ECF6","radialaxis":{"gridcolor":"white","ticks":"","linecolor":"white"}},"autotypenumbers":"strict","font":{"color":"#2a3f5f"},"ternary":{"baxis":{"gridcolor":"white","ticks":"","linecolor":"white"},"bgcolor":"#E5ECF6","caxis":{"gridcolor":"white","ticks":"","linecolor":"white"},"aaxis":{"gridcolor":"white","ticks":"","linecolor":"white"}},"annotationdefaults":{"arrowhead":0,"arrowwidth":1,"arrowcolor":"#2a3f5f"},"plot_bgcolor":"#E5ECF6","title":{"x":0.05},"scene":{"xaxis":{"gridcolor":"white","gridwidth":2,"backgroundcolor":"#E5ECF6","ticks":"","showbackground":true,"zerolinecolor":"white","linecolor":"white"},"zaxis":{"gridcolor":"white","gridwidth":2,"backgroundcolor":"#E5ECF6","ticks":"","showbackground":true,"zerolinecolor":"white","linecolor":"white"},"yaxis":{"gridcolor":"white","gridwidth":2,"backgroundcolor":"#E5ECF6","ticks":"","showbackground":true,"zerolinecolor":"white","linecolor":"white"}},"colorway":["#636efa","#EF553B","#00cc96","#ab63fa","#FFA15A","#19d3f3","#FF6692","#B6E880","#FF97FF","#FECB52"]},"data":{"barpolar":[{"type":"barpolar","marker":{"line":{"color":"#E5ECF6","width":0.5}}}],"carpet":[{"aaxis":{"gridcolor":"white","endlinecolor":"#2a3f5f","minorgridcolor":"white","startlinecolor":"#2a3f5f","linecolor":"white"},"type":"carpet","baxis":{"gridcolor":"white","endlinecolor":"#2a3f5f","minorgridcolor":"white","startlinecolor":"#2a3f5f","linecolor":"white"}}],"scatterpolar":[{"type":"scatterpolar","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}],"parcoords":[{"line":{"colorbar":{"ticks":"","outlinewidth":0}},"type":"parcoords"}],"scatter":[{"type":"scatter","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}],"histogram2dcontour":[{"colorbar":{"ticks":"","outlinewidth":0},"type":"histogram2dcontour","colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"contour":[{"colorbar":{"ticks":"","outlinewidth":0},"type":"contour","colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"scattercarpet":[{"type":"scattercarpet","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}],"mesh3d":[{"colorbar":{"ticks":"","outlinewidth":0},"type":"mesh3d"}],"surface":[{"colorbar":{"ticks":"","outlinewidth":0},"type":"surface","colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"scattermapbox":[{"type":"scattermapbox","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}],"scattergeo":[{"type":"scattergeo","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}],"histogram":[{"type":"histogram","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}],"pie":[{"type":"pie","automargin":true}],"choropleth":[{"colorbar":{"ticks":"","outlinewidth":0},"type":"choropleth"}],"heatmapgl":[{"colorbar":{"ticks":"","outlinewidth":0},"type":"heatmapgl","colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"bar":[{"type":"bar","error_y":{"color":"#2a3f5f"},"error_x":{"color":"#2a3f5f"},"marker":{"line":{"color":"#E5ECF6","width":0.5}}}],"heatmap":[{"colorbar":{"ticks":"","outlinewidth":0},"type":"heatmap","colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"contourcarpet":[{"colorbar":{"ticks":"","outlinewidth":0},"type":"contourcarpet"}],"table":[{"type":"table","header":{"line":{"color":"white"},"fill":{"color":"#C8D4E3"}},"cells":{"line":{"color":"white"},"fill":{"color":"#EBF0F8"}}}],"scatter3d":[{"line":{"colorbar":{"ticks":"","outlinewidth":0}},"type":"scatter3d","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}],"scattergl":[{"type":"scattergl","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}],"histogram2d":[{"colorbar":{"ticks":"","outlinewidth":0},"type":"histogram2d","colorscale":[[0.0,"#0d0887"],[0.1111111111111111,"#46039f"],[0.2222222222222222,"#7201a8"],[0.3333333333333333,"#9c179e"],[0.4444444444444444,"#bd3786"],[0.5555555555555556,"#d8576b"],[0.6666666666666666,"#ed7953"],[0.7777777777777778,"#fb9f3a"],[0.8888888888888888,"#fdca26"],[1.0,"#f0f921"]]}],"scatterternary":[{"type":"scatterternary","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}],"scatterpolargl":[{"type":"scatterpolargl","marker":{"colorbar":{"ticks":"","outlinewidth":0}}}]}},"showarrow":false,"margin":{"l":50,"b":50,"r":50,"t":60},"plot_bgcolor":"#f1f3f7","yaxis":{"showgrid":false,"showticklabels":false,"zeroline":false},"dragmode":"select"},
        {"editable":false,"responsive":true,"staticPlot":false,"scrollZoom":true},
    )
}

        });
    </script>
</div>

</body>
</html>
```


After the grid has been created it can be passed to the function `GetSystem()`, which then returns the system matrices for the state space representation in the continous time domain.


```julia
A, B, C, D = GetSystem(S2_L1);
```

We convert the matrices into the discrete time domain and create a discrete StateSpace object with the help of `ControlSystems`. This object can then represent the dynamics of the system for a given time interval using the function `lsim()`.


```julia
ts = 1e-5
Ad = exp(A*ts)
Bd = A \ (Ad - C) * B
sys_d = StateSpace(Ad, Bd, C, D, ts);
```


    UndefVarError: StateSpace not defined

    

    Stacktrace:

     [1] top-level scope

       @ In[8]:4


To use `lsim()` you need defined initial states `x0`, a time vector `t` and a input signal `u`. In our case we apply a jump to 250 V to the system.


```julia
ns = length(A[1,:]) # get num of states
ni = length(B[1,:]) # get num of inputs
t = collect(0:ts:0.1)
x0 = [0.0 for i = 1:ns]
u = [250.0 for i = 1:length(t)]
uu = [u for i = 1:ni ]
uuu = mapreduce(permutedims, vcat, uu);
```

`lsim()` now solves the difference equations for the given time steps and we can observe how the states evolve.


```julia
xout, _, _, _ = lsim(sys_d,uuu,t,x0=x0);
```


    UndefVarError: lsim not defined

    

    Stacktrace:

     [1] top-level scope

       @ In[10]:1


Here we plot the voltage across the capacitor in the first source.


```julia
layout = Layout(xaxis_title="Time in µs", yaxis_title="v_C / V")
p = plot(t, xout[2,:], layout)
```


    UndefVarError: Layout not defined

    

    Stacktrace:

     [1] top-level scope

       @ In[11]:1


### Random Grid Generation 

The manual creation of a grid can be time consuming because all parameters have to be specified. To quickly create different topologies with little effort, automatic generation can be used. t only needs the number of sources and loads, all other parameters are set by default or are randomly selected. An important feature is the creation of random node structures, where the parameters of the elements are chosen randomly.


```julia
S2_L2_FC = NodeConstructor(num_sources=2, num_loads=2);
```

    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit https://github.com/coin-or/Ipopt
    ******************************************************************************
    


To create the CM matrix randomly, the parameters `S2S_p`, `L2L_p` and `S2L_p` can also be used to set the connection probabilities. These indicate the probability with which a source/load will be connected to any other source/load. A fully connected grid would be obtained if all parameters are set to 1.


```julia
S2_L2_FC = NodeConstructor(num_sources=2, num_loads=2, S2S_p=1, S2L_p=1, L2L_p=1);
```

Let's check the CM matrix.


```julia
S2_L2_FC.CM
```




    4×4 Matrix{Float64}:
      0.0   1.0   2.0  3.0
     -1.0   0.0   4.0  5.0
     -2.0  -4.0   0.0  6.0
     -3.0  -5.0  -6.0  0.0



With a look into the parameter dict we also see that the parameters of the individual elements were randomly generated. The current policy for the sources is that an LC filter is always taken and the other filter types are chosen randomly.


```julia
S2_L2_FC.parameters["source"]
```




    2-element Vector{Any}:
     Dict{Any, Any}("L1" => 683065.1506262049, "C" => 770.6152043421689, "mode" => "Synchronverter", "fltr" => "LCL", "pwr" => 10000.0, "source_type" => "ideal", "L2" => 93413.60873776274, "R_C" => 9.924097510713226, "std_asy" => 2500.0, "σ" => 0.0…)
     Dict{Any, Any}("L1" => 455376.76708413655, "C" => 1155.9228065132534, "mode" => "Synchronverter", "fltr" => "LC", "pwr" => 15000.0, "source_type" => "ideal", "R_C" => 6.616065007142151, "std_asy" => 3750.0, "σ" => 0.0, "i_limit" => 52.183509595574726…)



For larger networks, of course, only the number of sources and loads can be handed over, so that the network structures are created on the basis of the default values. An important point here is that it is ensured that no subnets are created. By default, it is ensured that each element of the network has at least one connection to the other components of the network, so that no subnetworks are created.

Internally, this is done by checking for connections for each element. If these are not present, they are automatically created. For smaller networks it is advisable to specify a CM matrix, because otherwise usually too many connections are made than necessary. However, this is no longer noticeable with more than 10+ elements.


```julia
S5_L15 = NodeConstructor(num_sources=5, num_loads=15);
```

## Three-phase simulation

Until now, only single-phase grids have been created with NodeConstructor. However, the default value for the number of phases is 3, so we will now also consider the three-phase variant. For the simulation of the grid in 3 phases the system matrix is extended:

```math
\begin{equation}
A_{all} = \begin{pmatrix}
        A_{a} & 0 & 0 \\
        0 & A_{b} & 0 \\
        0 & 0 & A_{c} 
    \end{pmatrix}
\end{equation}
```

Lets see how this looks like in Julia.


```julia
S1_L2 = NodeConstructor(num_sources=1, num_loads=2);
```


```julia
S1_L2.parameters["grid"]["phase"]
```




    3




```julia
A, B, C, D = GetSystem(S1_L2)
ts = 1e-4
Ad = exp(A*ts)
Bd = A \ (Ad - C) * B
sys_d = StateSpace(Ad, Bd, C, D, ts);
```


    UndefVarError: StateSpace not defined

    

    Stacktrace:

     [1] top-level scope

       @ In[19]:5


We then collect a few useful variables and set the time horizon for the simulation.


```julia
ns = S1_L2.num_spp  # get num of states per phase
ni = S1_L2.num_sources # get num of inputs per phase
t = collect(0:ts:1);
```

Next we want to generate the three-phase input signals and repeat it for the number of sources:


```julia
# Stepfunction
u = sqrt(2)*[230, 0, -230]
uu = repeat(u, inner=ni) .* ones(length(t))';
```


```julia
# Sin wave
u = [230 * sin.(2*pi*t .- 2/3*pi*(i-1)) for i = 1:3]
uu = transpose(hcat(repeat(u[1], inner=[1,ni]),repeat(u[2], inner=[1,ni]),repeat(u[3], inner=[1,ni])));
```

Lets have a look:


```julia
layout = Layout(xaxis_title="Time in µs", yaxis_title="U in V")
input = 1

phase_a = scatter(x=t, y=uu[input+ni*0,:], mode="lines", name="Phase A")
phase_b = scatter(x=t, y=uu[input+ni*1,:], mode="lines", name="Phase B")
phase_c = scatter(x=t, y=uu[input+ni*2,:], mode="lines", name="Phase C")

plot([phase_a, phase_b, phase_c], layout)
```


    UndefVarError: Layout not defined

    

    Stacktrace:

     [1] top-level scope

       @ In[23]:1


Again, the discretized matrices can now be used to model the grid.


```julia
x0 = [0.0 for i = 1:ns*3]
xout, _, _, _ = lsim(sys_d,uu,t,x0=x0);
```


    UndefVarError: lsim not defined

    

    Stacktrace:

     [1] top-level scope

       @ In[24]:2


Now a state can be selected and the corresponding trajectories can be plotted.


```julia
state_list = GetStateIds(S1_L2)
```




    27-element Vector{String}:
     "source1_i_L1_a"
     "source1_v_C_filt_a"
     "source1_i_L2_a"
     "source1_v_C_cables_a"
     "cable1_i_L_a"
     "cable2_i_L_a"
     "load1_v_C_total_a"
     "load1_i_L_a"
     "load2_v_C_total_a"
     "source1_i_L1_b"
     "source1_v_C_filt_b"
     "source1_i_L2_b"
     "source1_v_C_cables_b"
     ⋮
     "load1_v_C_total_b"
     "load1_i_L_b"
     "load2_v_C_total_b"
     "source1_i_L1_c"
     "source1_v_C_filt_c"
     "source1_i_L2_c"
     "source1_v_C_cables_c"
     "cable1_i_L_c"
     "cable2_i_L_c"
     "load1_v_C_total_c"
     "load1_i_L_c"
     "load2_v_C_total_c"




```julia
state = 3
state_list = GetStateIds(S1_L2)

layout = Layout(xaxis_title="Time in µs", yaxis_title="$(state_list[state]) in V")

phase_a = scatter(x=t, y=xout[state+ns*0,:], mode="lines", name="Phase A")
phase_b = scatter(x=t, y=xout[state+ns*1,:], mode="lines", name="Phase B")
phase_c = scatter(x=t, y=xout[state+ns*2,:], mode="lines", name="Phase C")

plot([phase_a, phase_b, phase_c], layout)

```


    UndefVarError: Layout not defined

    

    Stacktrace:

     [1] top-level scope

       @ In[26]:4


### Access to the different states

A way to get the different states of our NodeConstructor is to use the function `GetStateIds()`.


```julia
state_list = GetStateIds(S1_L2)
```




    27-element Vector{String}:
     "source1_i_L1_a"
     "source1_v_C_filt_a"
     "source1_i_L2_a"
     "source1_v_C_cables_a"
     "cable1_i_L_a"
     "cable2_i_L_a"
     "load1_v_C_total_a"
     "load1_i_L_a"
     "load2_v_C_total_a"
     "source1_i_L1_b"
     "source1_v_C_filt_b"
     "source1_i_L2_b"
     "source1_v_C_cables_b"
     ⋮
     "load1_v_C_total_b"
     "load1_i_L_b"
     "load2_v_C_total_b"
     "source1_i_L1_c"
     "source1_v_C_filt_c"
     "source1_i_L2_c"
     "source1_v_C_cables_c"
     "cable1_i_L_c"
     "cable2_i_L_c"
     "load1_v_C_total_c"
     "load1_i_L_c"
     "load2_v_C_total_c"



The IDs created here are unique and can be used to access particular states. When creating the IDs, the sources are checked first in the order LCL, LC and then L. Then the cables are listed, which are also arranged in order. For the loads the order is RLC, LC, RL, L, RC, C and then R.

For the three-phase case, the state IDs are repeated and the respective phase is added.

These can then be accessed as follows:


```julia
state = 3
println(state_list[state+ns*0])
println(state_list[state+ns*1])
println(state_list[state+ns*2])
```

    source1_i_L2_a
    source1_i_L2_b
    source1_i_L2_c


Or:


```julia
state = "source1_i_L1"
idx_of_state = findall(x->occursin(state, x), state_list)
idx = idx_of_state
```




    3-element Vector{Int64}:
      1
     10
     19




```julia
state = "source1_i_L1_a"
idx_of_state = findall(x->occursin(state, x), state_list)
idx = idx_of_state
```




    1-element Vector{Int64}:
     1



The actions in the grid are also assigned unique IDs, here the sources are sorted in order. The IDs are output via the function `GetActionIds()`.


```julia
GetActionIds(S1_L2)
```




    3-element Vector{String}:
     "source1_u_a"
     "source1_u_b"
     "source1_u_c"




```julia

```
