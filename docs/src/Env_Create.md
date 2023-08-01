# How to build an environment using ElectricGrid

In the following the functionality of the `ElectricGridEnv` which creates an environment for simulation of power electronics-driven microgrids with ad-hoc data generation is described. 
It can be used to train and test reinforcement learing agents (e.g., from [juliareinforcementlearning.jl](https://juliareinforcementlearning.org/)), in particular for autonomous grid operation.
These agents can learn to handle different control tasks and can be compared to classical control approaches. 
The following topics are covered:

- ### Create an enviroment, 
- ### Setup an RL-env, 
- ### Detailed parameterisation.

The interactive content related to the section described here can be found in the form of a notebook [here](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/Env_Create_DEMO.ipynb).

## Create an Enviroment 

To create an environment consisting of an electrical power grid with two sources `num_sources = 2` (which could be for example an inverter fed by a PV plant or a wind turbine or battery) - supplying one load `num_sources = 1` two lines of code are sufficient like will be shown in the following.
An easy example is shown in the figure below, where a load (household) is supplied by 2 sources (inverters, fed PV plant and wind turbine) via two cables.

![](./assets/ExampleGrid.png)

A partially visible electrical circuit is displayed in the background as single phase diagram. 
By default a three-phase four wire system is created, for better visuablilty only one phase is shown here.
If it is not defined during the initialization of the env, all parameters (connections between the different sources and loads, parameters of the electric components,...) are drawn randomly, while a few are set to fixed values per default.
One of the latter would be for example the stepsize `ts`. 
The most important parameters will be investigated in the following.
A full list of all parameters which can be used for parametrisation can be found here (LINK TO PARAMETER DOCS).
After the initialization a step-wise interaction with the environment is possible.


To use the ElectricGrid.jl tool, the ElectricGrid package has to be loaded:


```julia
using ElectricGrid

env = ElectricGridEnv(num_sources = 2, num_loads = 1)
```
The created env  is conform to the interface like described in [juliareinforcementlearning.jl/env](https://juliareinforcementlearning.org/docs/How_to_write_a_customized_environment/).
Analog to this, the most important properties of the environment are:

 - is the environment terminated (`done`) which can happen in response to a limit violation or by reaching the number of interactions manually defined in the experiment,
 - the range of the inputs (`Action Space`) the environment accepts,
 - the range of the state (`State Space`) the environment returns,
 - the `Current State` of the environment.


In the following, the most important parts of the `ElectricGridEnv` will be investigated in detail. For this purpose, we will take a closer look at the generated `env`.

For better readability, the first five of the `Current State` of the environment are shown:


```julia
env.state[1:5]
```


    5-element Vector{Float64}:
     0.0
     0.0
     0.0
     0.0
     0.0


Since we have not yet interacted with the env, the `Current State` is still the initial - which by default is zero.
To learn how to interact with the environment applying different actions, see [Interaction with the Environment](https://upb-lea.github.io/ElectricGrid.jl/dev/Env_Interaction/).


If the state is not zero, the reset method can be used which sets the state to the internally defined `x0` (which consists of zeros per default).


The `Action Space` defines how many actions can be applied to the environment.
Since per default the env produces a three-phase system in the example shown above one action per phase per source is needed. 
This leads to 6 actions in total (2 sources with 3 phases each). 

The states represent physical values from the energy grid (for more detailed information see [The Nodecontructor - Theroy](https://upb-lea.github.io/ElectricGrid.jl/dev/NodeConstructor_Theory/)) 
The environment stores a specific index (`string`) per state.
That index highlights the physical context. 
Later it will be used during the control procedure (for more detailed information see [Mutlicontroller](https://upb-lea.github.io/ElectricGrid.jl/dev/RL_Classical_Controllers_Merge/)).


For the upper example the first five state indices for example are:
```julia
env.state_ids[1:5] # First 5 states
```


    5-element Vector{String}:
     "source1_i_L1_a"
     "source1_v_C_filt_a"
     "source1_i_L2_a"
     "source1_v_C_cables_a"
     "source2_i_L1_a"


The labels define to which source the state belongs and what it is. 
For example the first state is called `"source1_i_L1_a"`. That tells, it belongs to the first source (in the picture above the PV plant) and represents the current `i` through the incductor `L1` of phase `a`.
For example, this information can be used to control the current through the filter inductance (or to learn this control task).
All states of phase `a` of the example in the figure above (PV, wind turbine and houshold) are labeled in the circuit below (currents in red, voltages in green).

![](./assets/S2_L1_states.png)

The `action_ids` follow a similar concept:

```julia
env.action_ids[1:4] # First 4 actions
```


    4-element Vector{String}:
     "source1_u_a"
     "source2_u_a"
     "source1_u_b"
     "source2_u_b"


For example, `"source1_u_a"` belongs to the first source and represents the action `u` applied to phase `a`.

The `State Space` of the env typically from -1.0..1.0.
That means, the current `"source1_i_L1_a"` through the filter inductor in the example is normalized by the maximal current allowed to flow through the inductor.
If this parameter is not pre-defined it is set per default based on the filter layout happening in the env.
All (technical) parameters needed for the simulation are defined in the parameter dict (for detailed information see below and Default_Parameters.ipynb).
It can be investigated by:


```julia
env.nc.parameters
```


    Dict{Any, Any} with 4 entries:
      "source" => Any[Dict{Any, Any}("L1"=>0.0136613, "C"=>3.85308e-6, "mode"=>"Syn…
      "grid"   => Dict{Any, Any}("f_grid"=>50, "Δfmax"=>0.5, "fs"=>10000.0, "proces…
      "load"   => Any[Dict{Any, Any}("Z"=>24.0346-14.1138im, "C"=>0.000102027, "i_l…
      "cable"  => Any[Dict{Any, Any}("Cb"=>4.0e-7, "Lb"=>0.000264, "Rb"=>0.722, "C"…


The mentioned limit of the filter inductor current can be found using:


```julia
env.nc.parameters["source"][1]["i_limit"]
```


    16.18093320793015


Which returns the current limit (belonging to the inductor) of source one.
The voltage limit for normalization is depending on the filter capacitor and can be found using `env.nc.parameters["source"][1]["v_limit"]`.
The same concept holds for the cables and loads (parametrisation can be found in the parameter dict, too).

Since the `Action Space` is defined in a range -1.0..1.0, the actions are "normalized" by the DC-link voltage of the specific source. 
In the simulation the chosen action is multiplied by half of the DC-link voltage (and can be interpreted as modulation index in an electrical engineering context).
The DC-link voltage can be found in (or set via) the parameter dict, too:


```julia
env.nc.parameters["source"][1]["vdc"]
```


    800


The following figure shows the physical parameters of our example grid, which are also included in the parameter dict.

![](./assets/S2_L1_parameters.png)

## Setting Parameters

The wanted setting of the simulation can be adjusted using the parameter dict. As shown above, there are 4 keys: `"grid"`, `"source"`, `"cable"` and `"load"`. The most important parameters will be explained in the following:


### Grid
Defines the basic setting of the whole electrical power grid. The most important parameters are:

- `"f_grid"`: grid frequency (frequency of the voltage sinusoidal signal) (default: 50 Hz)
- `"phases"`: number of phases in the electric power grid (exclusive neutral). (default: 3)
- `"fs"`: sample frequency of the simulation. Every step the environment if simulated `ts = 1/fs` s forward (default: 10 kHz)
- `"v_rms"`: root mean square value of the basic grid voltage (default: 230 V)
- ...



```julia
env.nc.parameters["grid"]
```


    Dict{Any, Any} with 9 entries:
      "f_grid"        => 50
      "Δfmax"         => 0.5
      "fs"            => 10000.0
      "process_start" => 0.04
      "phase"         => 3
      "ramp_end"      => 0.04
      "ΔEmax"         => 5
      "pwr"           => 35000.0
      "v_rms"         => 230


### Source
Here, a list with all sources of the network can be found. Each source is in turn stored as a dict.


```julia
env.nc.parameters["source"]
```


    2-element Vector{Any}:
     Dict{Any, Any}("L1" => 0.013661303012524098, "C" => 3.853076021710845e-6, "mode" => "Synchronverter", "fltr" => "LCL", "pwr" => 5000.0, "source_type" => "ideal", "L2" => 0.0018682721747552546, "R_C" => 19.848195021426452, "std_asy" => 1250.0, "σ" => 0.0…)
     Dict{Any, Any}("L1" => 0.002276883835420683, "C" => 2.3118456130265067e-5, "mode" => "Synchronverter", "fltr" => "LC", "pwr" => 30000.0, "source_type" => "ideal", "R_C" => 3.3080325035710754, "std_asy" => 7500.0, "σ" => 0.0, "i_limit" => 104.36701919114945…)


A selection of the most important parameters that can be defined and configured via the parameter dict is given in the following:

- `"pwr"`: maximum apparent power of the source (default: random)
- `"source_type"`: type of the electric component sitting on the DC side, e.g., ideal (constant `vdc`), PV,... (default: ideal)
- `"control_type"`: defines whether the source is classically controlled or RL controlled (default = "classic")
- `"mode"`: allows to specify which control mode the source (default = "Droop", and classic control_type) (for more information, see Classic_Controllers_Notebook.iypnb)
- `"vdc"`: DC-link voltage, fixed if `"source_type"` is ideal (drawn random $U$[690, 800] V), otherwise drawn from function
- `"fltr"`: Filter type [L, LC, LCL] (default: random)
- `"i_limit"`: maximal allowed current flowing through the inductor(s) (default: calculated based on filter parameters)
- `"v_limit"`: maximal allowed voltage across the capacitor (default: calculated based on filter parameters).
- ...

The parameters of the first source are for example:


```julia
env.nc.parameters["source"][1]
```


    Dict{Any, Any} with 30 entries:
      "L1"           => 0.0136613
      "C"            => 3.85308e-6
      "mode"         => "Synchronverter"
      "fltr"         => "LCL"
      "pwr"          => 5000.0
      "source_type"  => "ideal"
      "L2"           => 0.00186827
      "R_C"          => 19.8482
      "std_asy"      => 1250.0
      "σ"            => 0.0
      "i_limit"      => 16.1809
      "v_rip"        => 0.01537
      "v_δ_set"      => 0.0
      "vdc"          => 800
      "τv"           => 0.002
      "k"            => 0
      "control_type" => "classic"
      "v_pu_set"     => 1.0
      "τf"           => 0.002
      ⋮              => ⋮


### Load
Here, a list with all loads of the network can be found. Each load is in turn stored as a dict.


```julia
env.nc.parameters["load"]
```


    1-element Vector{Any}:
     Dict{Any, Any}("Z" => 24.03463015163016 - 14.113790073837922im, "C" => 0.0001020274336341303, "i_limit" => 10.502997438247528, "L" => 0.22924428595808855, "R" => 32.322632462947404, "pwr" => 5693.8345453834845, "pf" => -0.862313846520903, "impedance" => "RLC")


Defines passive loads which can be pluged/connected to the grid. All combinations of the passive components of a resistor, inductor and capacitor ("impedance") can be chosen:

- `"impedance"`: type/circuit of the load (can be something out of [RLC, RL, RC, LC, R, L, C]). Components will be connected in parallel (default: random)
- `"pwr"`: apparant power drawn from the load assuming fixed vrms grid voltage and frequency and the number of phases (see parameters of "grid")
- `"pf"`: power factor resulting from chosen parameters
- `"Z"`: impedance resulting from chosen parameters
- `"R"`: load resistance
- ...


The parameters of the load, since only one is defined, are for example:


```julia
env.nc.parameters["load"][1]
```


    Dict{Any, Any} with 8 entries:
      "Z"         => 24.0346-14.1138im
      "C"         => 0.000102027
      "i_limit"   => 10.503
      "L"         => 0.229244
      "R"         => 32.3226
      "pwr"       => 5693.83
      "pf"        => -0.862314
      "impedance" => "RLC"


### Cable
Here, a list with all cables of the network can be found. Each cable is in turn stored as a dict.


```julia
env.nc.parameters["cable"]
```


    2-element Vector{Any}:
     Dict{Any, Any}("Cb" => 4.0e-7, "Lb" => 0.000264, "Rb" => 0.722, "C" => 0.0004, "i_limit" => 1.0e13, "v_limit" => 374.05948724768365, "len" => 1.0, "L" => 0.00025, "R" => 0.208)
     Dict{Any, Any}("Cb" => 4.0e-7, "Lb" => 0.000264, "Rb" => 0.722, "C" => 0.0004, "i_limit" => 1.0e13, "v_limit" => 374.05948724768365, "len" => 1.0, "L" => 0.00025, "R" => 0.208)


The cables are modelled using PI-models like shown in the figure above.
If no parameters are defined they are automatically designed based on the power flowing through the cable.


- `"len"`: length of the cable (default: random [1m, 1km])
- `"i_limit"`: maximal allowed current flowing through the inductor (default: calculated based on power flow and inductance)
- `"v_limit"`: since the capacitance is added to the one (if) defined in the source, the limit for the voltage can be found in the source parameter dict
- `"Cb"`: cable capacity coatings (default: 0.4 µF/km)
- `"Lb"`: cable operating inductor (default: 0.264 mH/km)
- `"Rb"`: cable AC resistor (default: 0.722 Ohm/km)
- `"C"`: cable capacity (default: `"len"`*`"Cb"`)
- `"L"`: cable inductance (default: `"len"`*`"Lb"`)
- `"R"`: cable resistance (default: `"len"`*`"Rb"`)

The parameters of the first cable are for example:


```julia
env.nc.parameters["cable"][1]
```


    Dict{Any, Any} with 9 entries:
      "Cb"      => 4.0e-7
      "Lb"      => 0.000264
      "Rb"      => 0.722
      "C"       => 0.0004
      "i_limit" => 1.0e13
      "v_limit" => 374.059
      "len"     => 1.0
      "L"       => 0.00025
      "R"       => 0.208


For all information about the parameters which can be defined investigate the parameter dict itself or [have a look into the documentation](https://upb-lea.github.io/ElectricGrid.jl/dev/Default_Parameters/).

The next step would be to interact with the env to run an experiment.
To learn how to do this, see [Interaction with the Environment](https://upb-lea.github.io/ElectricGrid.jl/dev/Env_Interaction/).
