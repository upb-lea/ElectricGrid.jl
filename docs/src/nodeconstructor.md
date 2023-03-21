# NodeConstructor Demo

The NodeConstructor has the goal to convert the physical system, here the grid, into a mathematical model. This mathematical model can then be used to simulate the behavior of the network. This notebook addresses the following points:

- ### Introduction to the NodeConstructor
- ### Representation of the physical grid
- ### Building of the state-space matrices


## Introduction

![Illustration of where the NodeConstructor is located](./assets/OverviewJEG.png)

The NodeConstructor is part of the enviroment (see graphic). Its function is to create the mathematical model of the grid to be simulated. The inputs for the NodeConstructor are the specifications of the grid, which will be discussed in more detail subsequently. Within the NodeConstructor, an ordinary differential equation (ODE) system is created based on the underlying physical models and properties of the individual components. Based on this ODE system, the necessary matrices can be extracted.

In the following, we will discuss how the mathematical model of an electrical power grid can be generated using the NodeConstructor.

## Modelling a physical grid

The following graphic shows a simple example grid:


![Illustration of a simple network and its electrical equivalent circuit diagram.](./assets/ExampleGrid.png)

Two sources (wind turbine on the left and PV arrays on the right), as well as a load (exemplary of a household which is composed of various objects) can be seen in this. In addition, the sources are connected to the load via cables. For the sources, we assume that the voltage is provided by a DC link, so we can assume that the voltage can be modeled by a voltage source and a filter. Loads are divided into active loads, which, based on an internal function, can take power from the network, and passive loads, which are described by their resistive, capacitive and/or inductive components and thus dissipate a certain power. Cable modeling relies on the Pi equivalent circuit like depicted in the fiure below.

If you abstract the above figure and insert the respective components according to the description, you get the following single-phase circuit:

![Abstraction of the previously presented example grid.](./assets/ExampleGridCurcuit.png)

With the help of this example, we will now discuss the derivation of a mathematical model and what information the NodeConstructor will later need to generate this model automatically.

### Fundamentals of electrical engineering

If we annotate the above equivalent circuit, we get the following representation:

![Equivalent curcuit with annotations.](./assets/ExampleGridCurcuit_Annotations_Kirchhof.png)

In this, the components are annotated and the currents are noted with the technical flow direction. For the sake of clarity, the voltages are only indicated by green arrows, but have the same name as the corresponding component to which they are applied.

Based on the components within the equivalent circuit diagram, we can now construct the differential equations. Togehter with Kirchhoff's laws

$$
\begin{align}
    \sum_k^n i_k(t) &= 0,\\
    \sum_k^n u_n(t) &= 0
\end{align}
$$

and the component equations

$$
\begin{align}
    u(t) &= R i(t),\\
    i(t) &= C\dot u(t),\\
    u(t) &= L\dot i(t).
\end{align}
$$

This gives the following equation system:

$$
\begin{align}
    \text{M1:}\;& 0 = u_{in1}  - u_{R11} - u_{L11} - u_{C11} - u_{RC1}, \tag{1}\label{equ:1}\\
    \text{M2:}\;& 0 = u_{C11} + u_{RC1} - u_{Cb11},\tag{2}\label{equ:2}\\
    \text{M3:}\;& 0 = u_{Cb11} - u_{Rb1} - u_{Lb1} - u_{Cb21},\tag{3}\label{equ:3}\\
    \text{M4:}\;& 0 = u_{RL} - u_{LL}, \tag{4}\label{equ:4}\\
    \text{M5:}\;& 0 = u_{Cb22} + u_{Rb2} + u_{Lb2} - u_{Cb12}, \tag{5}\label{equ:5}\\
    \text{M6:}\;& 0 = u_{Cb12} + u_{R22} + u_{L22} - u_{C12} - u_{RC2}, \tag{6}\label{equ:6}\\
    \text{M7:}\;& 0 = u_{RC2} + u_{C12} + u_{L12} + u_{R12} - u_{in2}, \tag{7}\label{equ:7}\\
    \text{N1:}\;& 0 = i_{11} - i_{RC1} - i_{Cb11} - i_{Rb1}, \tag{8}\label{equ:8}\\
    \text{N2:}\;& 0 = i_{Rb1} - i_{Cb21} - i_{CL} - i_{RL} - i_{LL} - i_{Cb22} + i_{Rb2}, \tag{9}\label{equ:9}\\
    \text{N3:}\;& 0 = i_{22} - i_{Rb2} - i_{Cb21}, \tag{10}\label{equ:10}\\
    \text{N4:}\;& 0 = i_{12} - i_{RC2} - i_{22}. \tag{11}\label{equ:11}
\end{align}
$$


For the sake of clarity, the time dependencies are omitted here. For node 2 the connected capacitances are summed up, because there are several capacitors connected in parallel, e.g.:

\begin{equation}
    C_{SL} = C_{b1} + C_L + C_{b2}
\end{equation}

We then call this capacitor $C_{SL}$ through which the current $i_{CSL}$ flows and the voltage $u_{CSL}$ is applied. Now we insert the component equations into the equations \ref{equ:1} to \ref{equ:11} and rearrange them,  so that an ODE system is obtained:

\begin{align}
    \dot i_{11} &= \frac{1}{L_{11}} \left( u_{in1}  - R_{11} i_{11} - R_{RC1} i_{RC1} - u_{C12} \right), \\
    \dot u_{C11} &= \frac{1}{C_{11}R_{C11}} \left(u_{Cb11} - u_{RC1} \right), \\
    \dot i_{Rb1} &= \frac{1}{L_{b1}} \left( u_{Cb11} - i_{Rb1} R_{b1} \right),\\
    \dot i_{LL} &= \frac{1}{L_{L}} u_{CL}\\
    \dot i_{Rb2} &= \frac{1}{L_{b2}} \left( u_{Cb12} - i_{Rb2} R_{b2} - u_{Cb22}\right),\\
    \dot i_{22} &= \frac{1}{L_{22}} \left( u_{C12}  - i_{22} (R_{C2}+R_{22}) - R_{RC1} i_{RC2} - u_{Cb12} \right), \\
    \dot i_{12} &= \frac{1}{L_{12}} \left( u_{in2}  - i_{12} (R_{12} + R_{RC1}) - R_{RC1} i_{C1} - u_{C12} \right), \\
    \dot u_{Cb11} &= \frac{1}{C_{b1}} \left(i_{11} - i_{Rb1} - \frac{1}{R_{C1}} u_{Cb11} + \frac{1}{R_{C1}} u_{C11}\right)\\
    \dot u_{CSL} &= \frac{1}{C_{SL}} \left( i_{Rb1} + i_{Rb2} - \frac{1}{R_{L}} u_{CSL} - i_{LL} \right)\\
    \dot u_{Cb12} &= \frac{1}{C_{b2}} \left(i_{22} - i_{Rb2}\right).\\
    \dot u_{C12} &= \frac{1}{C_{12}} \left(i_{12} - i_{22}\right).\\
\end{align}

These differential equations describe the system completely and could be solved for exampled using numerical solvers like Forward Euler or Runge-Kutta. 

### Forming the system matrices from the ODEs

The differential equations found are now to be converted into the common state-space representation:

\begin{align}
    \dot{\vec{x}}(t) &= \mathbf{A} \vec{x}(t) + \mathbf{B} \vec{u}(t)\\
    \vec{y}(t) &= \mathbf{C} \vec{x}(t) + \mathbf{D} \vec{u}(t)
\end{align}

Since we are interested in how our grid evolves over time, the first equation is most relevant to our application. $x(t)$ describes in general the state vector, $\dot x(t)$ the changes of this state vector. The change of the system is described on the one hand by the system dynamics, which is expressed by the matrix $\mathbf{A}$ and on the other hand by the input signals $u(t)$ acting on the system combined with the matrix $\mathbf{B}$.

For our example, we can find the following state vector:

\begin{equation}
    \vec{x}(t) = 
    \begin{pmatrix}
        i_{11} & u_{C11} & u_{Cb11}  & i_{12} & u_{C12} & i_{22} & u_{Cb12} & i_{Rb1} & i_{Rb2} & u_{CSL} & i_{LL}
    \end{pmatrix}^\top
\end{equation}

and the following input vector:

\begin{equation}
    \vec{u}(t) = 
    \begin{pmatrix}
        u_{in1} & u_{in2}
    \end{pmatrix}^\top.
\end{equation}

It should be noted that the order of the states here is first sources, then cables and finally the loads. 

This leads for the investigated example to the following matrices:

\begin{equation}
    \mathbf{A} = 
    \begin{pmatrix}
        -\frac{R_{11}}{L_{11}} & 0 & -\frac{1}{L_{11}} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
        0 & -\frac{1}{C_{11} R_{C1}} & \frac{1}{C_{11} R_{C1}} & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
        \frac{1}{C_{b1}} & \frac{1}{C_{b1} R_{C1}} & -\frac{1}{C_{b1} R_{C1}} & 0 & 0 & 0 & 0 & -\frac{1}{C_{b1}} & 0 & 0 & 0\\
        0 & 0 & 0 & -\frac{R_{12} + R_{C2}}{L_{12}} & -\frac{1}{L_{12}} & \frac{R_{C2}}{L_{12}} & 0 & 0 & 0 & 0 & 0\\
        0 & 0 & 0 & \frac{1}{C_{12}} & 0 & -\frac{1}{C_{12}} & 0 & 0 & 0 & 0 & 0\\
        0 & 0 & 0 & \frac{R_{C2}}{L_{22}} & \frac{1}{L_{22}} & -\frac{R_{22} + R_{C2}}{L_{22}} & -\frac{1}{L_{22}} &0&0&0&0\\
        0 & 0 & 0 & 0 & 0 & \frac{1}{C_{b2}} & 0 & 0 & -\frac{1}{C_{b2}} & 0 & 0\\
        0 & 0 & \frac{1}{L_{b1}} & 0 & 0 & 0 & 0 & -\frac{R_{b1}}{L_{b1}} & 0 & -\frac{1}{L_{b1}} & 0\\
        0 & 0 & 0 & 0 & 0 & 0 & \frac{1}{L_{b2}} & 0 & -\frac{R_{b2}}{L_{b1}} & -\frac{1}{L_{b2}} & 0\\
        0 & 0 & 0 & 0 & 0 & 0 & 0 & \frac{1}{C_{SL}} & \frac{1}{C_{SL}} & -\frac{1}{R_{L} C_{SL}} & -\frac{1}{C_{SL}}\\
        0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & 0 & \frac{1}{L_{L}} & 0\\
    \end{pmatrix}
\end{equation}

and

\begin{equation}
    \mathbf{B} = 
    \begin{pmatrix}
        \frac{1}{L_{11}} & 0 \\
        0 & 0 \\
        0 & 0 \\
        0 & \frac{1}{L_{12}} \\
        0 & 0 \\
        0 & 0 \\
        0 & 0 \\
        0 & 0 \\
        0 & 0 \\
        0 & 0 \\
        0 & 0 \\
    \end{pmatrix}
\end{equation}

### Solving the ODE system

Theory of exact discritization

We now have the state-space matrices we need to simulate the grid.  There are several tools to solve the differential equations. We use the `ControlSystem` package ([Source](https://juliacontrol.github.io/ControlSystems.jl/latest/man/creating_systems/#ss-Creating-State-Space-Systems)), to be more precise the `lsim()` function. (AST Lecture 3, Slide A_d 47 folgende, Quelle)


```julia
using ControlSystemsBase
using LinearAlgebra
using PlotlyJS
```

Next we define some the values of the components:


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

Now we construct the matrices:


```julia
A = zeros((11,11))
A[1,1] = -R/L
A[1,3] = -1/L
A[2,2] = -1/(C*R_c)
A[2,3] = 1/(C*R_c)
A[3,1] = 1/C_b
A[3,2] = 1/(R_c*C_b)
A[3,3] = -1/(R_c*C_b)
A[3,8] = -1/C_b
A[4,4] = -(R+R_c)/L
A[4,5] = -1/L
A[4,6] = R_c/L
A[5,4] = 1/C
A[5,6] = -1/C
A[6,4] = R_c/L
A[6,5] = 1/L
A[6,6] = -(R+R_c)/L
A[6,7] = -1/L
A[7,6] = 1/C_b
A[7,9] = -1/C_b
A[8,3] = 1/L_b
A[8,8] = -R_b/L_b
A[8,10] = -1/L_b
A[9,7] = 1/L_b
A[9,9] = -R_b/L_b
A[9,10] = -1/L_b

C_SL = C_b + C_l + C_b

A[10,8] = 1/C_SL
A[10,9] = 1/C_SL
A[10,10] = -1/(R_l*C_SL)
A[10,11] = -1/C_SL
A[11,10] = 1/L_l;
```


```julia
B = zeros((11,2))

B[1,1] = 1/L
B[4,2] = 1/L;
```


```julia
C = Diagonal(ones(11));
```


```julia
D = 0;
```

We convert the matrices into the discrete time domain and create a discrete StateSpace object with the help of `ControlSystems`. This object can then represent the dynamics of the system for a given time interval using the function `lsim()`.


```julia
ts = 1e-5
Ad = exp(A*ts)
Bd = A \ (Ad - C) * B
sys_d = StateSpace(Ad, Bd, C, D, ts);
```


```julia
ns = length(A[1,:]) # get num of states
ni = length(B[1,:]) # get num of inputs
t = collect(0:ts:0.1)
x0 = [0.0 for i = 1:ns]
u = [230.0 for i = 1:length(t)]
uu = [u for i = 1:ni ]
uuu = mapreduce(permutedims, vcat, uu);
```

To use `lsim()` you need defined initial states `x0`, a time vector `t` and a input signal `u`. In our case we apply a jump to 230 V to the system.


```julia
xout, _, _, _ = lsim(sys_d,uuu,t,x0=x0);
```


```julia
layout = Layout(xaxis_title="Time in µs", yaxis_title="v_C / V")
p = plot(t, xout[5,:], layout)
```

## Automatic generation of grids

The creation of a grid always follows a certain pattern and can therefore be automated. In this way, you can avoid errors that may arise from the handwritten rearrangement of the equations. It is therefore obvious to write an automatism for the generation, which generates the matrices on the basis of input parameters.

In this part of the notebook the above example shall be reproduced with the help of the NodeConstructor.

An important point is the information about how the components are connected by which cable. To pass this information to the NodeConstructor we introduce the connection matrix (`CM`).

### Connection matrix

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

We now import the JuliaElectricGrid package:


```julia
# Import the package
using JEG
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

After the grid has been created it can be passed to the function `GetSystem()`, which then returns the matrices for the state space representation in the continous time domain.


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

Here we plot the voltage across the capacitor in the first source.


```julia
layout = Layout(xaxis_title="Time in µs", yaxis_title="v_C / V")
p = plot(t, xout[2,:], layout)
```


```julia

```


```julia

```


```julia

```


```julia

```

Manuell parameterization

.... Split up here .....

Automatic generation

## Automatic generation of the grids

An important feature is the creation of random node structures, where the parameters of the elements are chosen randomly. For fully connected structures, this can be generated, for example, using the parameters `S2S_p`, `L2L_p` and `S2L_p`. These indicate the connection probability of a source/load with any other source/load. If these parameters are set to 1, a fully connected node is generated.


```julia
S2_L2_FC = NodeConstructor(num_sources=2, num_loads=2, S2S_p=1, S2L_p=1);
```

Let's check the CM matrix.


```julia
S2_L2_FC.CM
```

With a look into the parameter dict we also see that the parameters of the individual elements were randomly generated. The current policy for the sources is that an LC filter is always taken and the other filter types are chosen randomly.


```julia
S2_L2_FC.parameters["source"]
```

For larger networks, of course, only the number of sources and loads can be handed over, so that the network structures are created on the basis of the default values. An important point here is that it is ensured that no subnets are created. By default, it is ensured that each element of the network has at least one connection to the other components of the network, so that no subnetworks are created.

Internally, this is done by checking for connections for each element. If these are not present, they are automatically created. For smaller networks it is advisable to specify a CM matrix, because otherwise usually too many connections are made than necessary. However, this is no longer noticeable with more than 10+ elements.


```julia
S5_L15 = NodeConstructor(num_sources=5, num_loads=15);
```

## Three-phase simulation

Until now, only single-phase grids have been created with NodeConstructor. However, the default value for the number of phases is 3, so we will now also consider the three-phase variant.


```julia
S5_L15 = NodeConstructor(num_sources=1, num_loads=2);
```


```julia
S5_L15.parameters["grid"]["phase"]
```


```julia
A, B, C, D = GetSystem(S5_L15)
ts = 1e-4
Ad = exp(A*ts)
Bd = A \ (Ad - C) * B
sys_d = StateSpace(Ad, Bd, C, D, ts);
```

We then collect a few useful variables and set the time horizon for the simulation.


```julia
ns = S5_L15.num_spp  # get num of states per phase
ni = S5_L15.num_sources # get num of inputs per phase
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

Again, the discretized matrices can now be used to model the grid.


```julia
x0 = [0.0 for i = 1:ns*3]
xout, _, _, _ = lsim(sys_d,uu,t,x0=x0);
```

Now a state can be selected and the corresponding trajectories can be plotted.


```julia
state_list = GetStateIds(S5_L15)
```


```julia
state = 3
state_list = GetStateIds(S5_L15)

layout = Layout(xaxis_title="Time in µs", yaxis_title="$(state_list[state]) in V")

phase_a = scatter(x=t, y=xout[state+ns*0,:], mode="lines", name="Phase A")
phase_b = scatter(x=t, y=xout[state+ns*1,:], mode="lines", name="Phase B")
phase_c = scatter(x=t, y=xout[state+ns*2,:], mode="lines", name="Phase C")

plot([phase_a, phase_b, phase_c], layout)

```

## Access to the different states

A way to get the different states of our NodeConstructor is to use the function `GetStateIds()`.


```julia
state_list = GetStateIds(S5_L15)
```

The IDs created here are unique and can be used to access particular states. When creating the IDs, the sources are checked first in the order LCL, LC and then L. Then the cables are listed, which are also arranged in order. For the loads the order is RLC, LC, RL, L, RC, C and then R.

For the three-phase case, the state IDs are repeated and the respective phase is added.

These can then be accessed as follows:


```julia
state = 3
println(state_list[state+ns*0])
println(state_list[state+ns*1])
println(state_list[state+ns*2])
```

Or:


```julia
state = "source1_i_L1"
idx_of_state = findall(x->occursin(state, x), state_list)
idx = idx_of_state
```


```julia
state = "source1_i_L1_a"
idx_of_state = findall(x->occursin(state, x), state_list)
idx = idx_of_state
```

The actions in the grid are also assigned unique IDs, here the sources are sorted in order. The IDs are output via the function `GetActionIds()`.


```julia
GetActionIds(S5_L15)
```


```julia

```
