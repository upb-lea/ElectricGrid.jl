# ElectricGridEnv Demo

This notebook is intended to show the functionality of the ElectricGridEnv which creates an environment for simulation of power-electronic driven microgrids with ad-hoc data generation.
It can be used to to train and test reinforcement learing agents (e.g., from 
https://juliareinforcementlearning.org/).
These agents can learn to handle different control tasks and can be compared to classical control approaches.

The dynmaic bahaviour of the envorinment is simulated using linear state-space systems.
It interacts step-wise with the agent/controller like shown in the figure below.
Based on the input/action `u` at timestep `k` the state `x` is calculated.

  
![](assets/RL_env.png)


To use the JEG tool the JEG package has to be loaded:


```julia
using JEG
```




## Simplest initialisation
The easiest way to initialize an environment is as follows:


```julia
env = ElectricGridEnv(num_sources = 2, num_loads = 1)
```

    
    ******************************************************************************
    This program contains Ipopt, a library for large-scale nonlinear optimization.
     Ipopt is released as open source code under the Eclipse Public License (EPL).
             For more information visit https://github.com/coin-or/Ipopt
    ******************************************************************************
    
    


# ElectricGridEnv

## Traits

| Trait Type        |                                            Value |
|:----------------- | ------------------------------------------------:|
| NumAgentStyle     |          ReinforcementLearningBase.SingleAgent() |
| DynamicStyle      |           ReinforcementLearningBase.Sequential() |
| InformationStyle  | ReinforcementLearningBase.ImperfectInformation() |
| ChanceStyle       |           ReinforcementLearningBase.Stochastic() |
| RewardStyle       |           ReinforcementLearningBase.StepReward() |
| UtilityStyle      |           ReinforcementLearningBase.GeneralSum() |
| ActionStyle       |     ReinforcementLearningBase.MinimalActionSet() |
| StateStyle        |     ReinforcementLearningBase.Observation{Any}() |
| DefaultStateStyle |     ReinforcementLearningBase.Observation{Any}() |

## Is Environment Terminated?

No

## State Space

`ReinforcementLearningBase.Space{Vector{IntervalSets.ClosedInterval{Float64}}}(IntervalSets.ClosedInterval{Float64}[-1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0])`

## Action Space

`ReinforcementLearningBase.Space{Vector{IntervalSets.ClosedInterval{Float64}}}(IntervalSets.ClosedInterval{Float64}[-1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0])`

## Current State

```
[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
```



This creates an environment consisting of an electrical power grid with two sources `num_sources = 2` (which could be for example an inverter fed by a PV plant or a wind turbine or battery) - supplying one load `num_sources = 1`.
An easy exemplary example is shown in the figure below, where a load (household) is supplied by 2 sources (inverters, fed PV plant and wind turbine) via two transmission lines.

![](assets/ExampleGrid.png)

For better visuablilty an exemplary shaded electircal circut in the background is displayed as single phase diagram.
(By default a three-phase four wire system is created).
If it is not defined during the initialization of the env, all parameters (connections between the different sources and loads, parameters of the electric components,...) are drawn randomly, while a few are set to fixed value per default.
One of the latter would be for example the stepsize `ts`. After the initialization a step-wise interaction with the environment is possible. 
As can be seen in the first picture, an action can be selected and the env can be executed with it. 
Based on that action `u_k` and the internal state-space system (defined depending on the electric components - for more information about the odernary differential equation,... see NodeConstructor_DEMO.ipynb) the system is evolved for one timestep and the new states `x_k+1` of the system are calulated.

### States and actions

First, the current states of the environment are checked:



```julia
env.state
```


    30-element Vector{Float64}:
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0
     ⋮
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0
     0.0


If the state is not zero, but should be in the beginning, the reset method can be used which sets the state to the internally defined `x0` (which consists of zeros per default).
If we do not want to not start from zero, we could set `x0` in the initialisation of the env:


```julia
using ReinforcementLearning
env.x0 = 0.1 * ones(length(env.state_space))
reset!(env)
env.state
```


    30-element Vector{Float64}:
     0.1
     0.1
     0.1
     0.1
     0.1
     0.1
     0.1
     0.1
     0.1
     0.1
     ⋮
     0.1
     0.1
     0.1
     0.1
     0.1
     0.1
     0.1
     0.1
     0.1


To interact with the env, it has to be figured out how many actions are needed. Therefore the length of the action space can be checked:





```julia
n_a = length(env.action_space)
```


    6


The six action requested by the environment belong to the 2 sources. Since per default the env produces a three-phase system we need one action per phase per source -> 6 actions.
To excite the env by an action the following command can be used:


```julia
env([0.2, 0.2, 0.2, 0.3, 0.3, 0.3])
```


    30-element Vector{Float64}:
     -0.16394393819718883
      0.008632257605152691
      2.0598878480930146e-8
     -0.1605061258652815
      0.027306112731665413
     -0.0016519374537311848
     -0.001854826977574919
      0.0002007141020688965
      0.031641836651051955
      0.13796560159583626
      ⋮
      0.00863225760515267
      2.0598878527674185e-8
     -0.1605061258652815
      0.027306112731665506
     -0.0016519374537312126
     -0.001854826977574919
      0.0002007141020688686
      0.031641836651051844
      0.13796560159583626


Here, the first source got an action of `0.2` to all three phases, while the second source got an action of `0.3` to all three phases.
As can be seen, the states have changed from 0.1 to different values.
To get a little bit more intuition about the different states, the state_ids can be investigated:


```julia
env.state_ids
```


    30-element Vector{String}:
     "source1_i_L1_a"
     "source1_v_C_filt_a"
     "source1_v_C_cables_a"
     "source2_i_L1_a"
     "source2_v_C_cables_a"
     "cable1_i_L_a"
     "cable2_i_L_a"
     "cable3_i_L_a"
     "load1_v_C_total_a"
     "load1_i_L_a"
     ⋮
     "source1_v_C_filt_c"
     "source1_v_C_cables_c"
     "source2_i_L1_c"
     "source2_v_C_cables_c"
     "cable1_i_L_c"
     "cable2_i_L_c"
     "cable3_i_L_c"
     "load1_v_C_total_c"
     "load1_i_L_c"


The labels define to which source the state belongs and what it is about. 
For example the first state is called `"source1_i_L1_a"`. That tells, it belongs to the first source (in the picture above the PV plant) and represents the current `i` through the incductor `L1` of phase `a`.
For example, this information can be used to control the current through the filter inductance (or to learn this control task).

An example for a single phase 3 Bus grid consisting of 2 sources (LC and LCL filter), 3 cables and 1 RLC load is shown in the following figure highligthing some of the defineable parameters and states.
This example only shows the equivalent electrical circuit diagramm but is similar to the example with the PV, windturbine and electric car. The only difference here is that the two sources are connected, too. 

![](assets/ExampleGrid3.png)

A few states are labeled in red, while a few parameters defined in the parameter dict are marked in green.


```julia
env.state_space
```


    Space{Vector{IntervalSets.ClosedInterval{Float64}}}(IntervalSets.ClosedInterval{Float64}[-1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0  …  -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0, -1.0..1.0])


Since the state space of the env tells, that it ranges from -1.0..1.0, the current `"source1_i_L1_a"` through the filter inductor in the example is normalized by the maximal current allowed to flow through the inductor.
If this parameter is not defined it it set per default based on the filter layout happening in the env.



All (technical) parameters needed for the simulation are defined in the parameter dict (for more detailed information see below and NodeConstructor_DEMO.ipynb).
It can be investigated by:




```julia
env.nc.parameters

```


    Dict{Any, Any} with 4 entries:
      "source" => Any[Dict{Any, Any}("L1"=>0.00170766, "C"=>3.08246e-5, "mode"=>"Sy…
      "grid"   => Dict{Any, Any}("f_grid"=>50, "Δfmax"=>0.005, "fs"=>10000.0, "proc…
      "load"   => Any[Dict{Any, Any}("Z"=>8.64561-6.22694im, "C"=>0.000255711, "L"=…
      "cable"  => Any[Dict{Any, Any}("Cb"=>4.0e-7, "Lb"=>0.000264, "Rb"=>0.722, "C"…


The limit of the filter inductor current can be found using:


```julia
env.nc.parameters["source"][1]["i_limit"]
```


    139.15602558819927


Which returns the current limit (belonging to the inductor) of source one.
The voltage limit for normalization is depending on the filter capacitor and calulated - if not defined - based on the capacitance.
The same concept holds for the cables and loads (parametrization can be found in the parameter dict, too).

Since the action space is defined in a range -1.0..1.0, the actions are "normalized" by the DC-link voltage of the specific source. 
In the simulation the chosen action is multiplied by half of the DC-link voltage (and can be interpreted as modulation index in an electrical engineering context).
The DC-link voltage can be found in (or set via) the parameter dict, too:


```julia
env.nc.parameters["source"][1]["vdc"]
```


    800


### Parameters

The wanted setting of the simulation can be adjusted using the parameter dict. 
The most important parts will be investigated in the following.
Like already shown above the parameter dict splits up into different parts describing different sections of the grid: gird, source, load, cable.

Now the most important parts of the parameter dict are investigated in more detail.

#### Grid
Defines the basic setting of the whole electrical power grid. The moste important parameters are:

- `"f_grid"`: grid frequency (frequency of the sine wave) (default: 50 Hz)
- `"phases"`: number of phases in the electric power grid (exclusive neutral). (default: 3)
- `"fs"`: sample frequency of the simulation. Every step the environment if simulated `1/fs` s forward (default: 10 kHz)
- `"v_rms"`: root mean square value of the basic grid voltage (default: 230 V)
- ...



```julia
env.nc.parameters["grid"]
```


    Dict{Any, Any} with 9 entries:
      "f_grid"        => 50
      "Δfmax"         => 0.005
      "fs"            => 10000.0
      "process_start" => 0.04
      "phase"         => 3
      "ramp_end"      => 0.04
      "ΔEmax"         => 0.05
      "pwr"           => 85000.0
      "v_rms"         => 230


#### Source
A selection of the most important parameters that can be defined and configured via the parameter dict is given in the following:

- `"pwr"`: maximum aparent power of the source (default: random)
- `"source_type"`: type of the electric component sitting on the DC side, e.g. ideal (constant `vdc`), PV,... (default: ideal)
- `"control_type"`: defines whether the source is classically controlled or RL controlled (default = "classic")
- `"mode"`: allows to specify which control mode the source (default = "Droop", and classic control_type) (for more information, see Classic_Controllers_Notebook.iypnb)
- `"vdc"`: DC-link voltage, fixed if `"source_type"` is ideal (drawn random $U$[690, 800] V), otherwise drawn from function
- `"fltr"`: Filter type [L, LC, LCL] (default: random)
- `"i_limit"`: maximal allowed current flowing through the inductor(s) (default: calculated based on filter parameters)
- `"v_limit"`: maximal allowed voltage across the capacitor (default: calculated based on filter parameters).
- ...




```julia
env.nc.parameters["source"][1]
```


    Dict{Any, Any} with 28 entries:
      "L1"           => 0.00170766
      "C"            => 3.08246e-5
      "mode"         => "Synchronverter"
      "fltr"         => "LC"
      "pwr"          => 40000.0
      "source_type"  => "ideal"
      "R_C"          => 2.48102
      "std_asy"      => 10000.0
      "σ"            => 0.0
      "i_limit"      => 139.156
      "v_rip"        => 0.01537
      "v_δ_set"      => 0.0
      "vdc"          => 800
      "τv"           => 0.002
      "k"            => 0
      "control_type" => "classic"
      "v_pu_set"     => 1.0
      "τf"           => 0.002
      "i_rip"        => 0.15
      ⋮              => ⋮


#### Load
Defines passive loads which can be pluged/connected to the grid. All combinations of the passive components of a resistor, inductor and capacitor ("impedance") can be chosen:

- `"impedance"`: type/circuit of the load (can be something out of [RLC, RL, RC, LC, R, L, C]). Components will be connected in parallel (default: random)
- `"pwr"`: apparant power of drawn from the load assuming fixed vrms grid voltage and frequency and the number of phases (see parameters of "grid")
- `"pf"`: power factor resulting from chosen parameters
- `"Z"`: impedance resulting from chosen parameters
- `"R"`: load resistance
- ...



```julia
env.nc.parameters["load"]
```


    1-element Vector{Any}:
     Dict{Any, Any}("Z" => 8.645614985778817 - 6.226943873467353im, "C" => 0.00025571060022957864, "L" => 0.12491902487968094, "R" => 13.130527865555832, "pwr" => 14894.91016754301, "pf" => -0.811440822661242, "impedance" => "RLC")


#### Cable
The cables are modelled using PI-models like shown in the figure above.
If no parameters are defined they are chosen based on the power flowing through the cable.


- `"len"`: length of the cable (default: random [1m, 1km])
- `"i_limit"`: mamimal allowed current flowing through the inductor (default: calculated based on power flow and inductance)
- `"v_limit"`: since the capacitance it added to the one (if) defined in the source, the limit for the voltage can be found in the source parameter dict
- `"Cb"`: cable capacity coatings (default: 0.4 µF)
- `"Lb"`: cable operating inductor (default: 0.264 mH)
- `"Rb"`: cable AC resistor (default: 0.722 Ohm)
- `"C"`: cable capacity (default: `"len"`*`"Cb"`)
- `"L"`: cable inductance (default: `"len"`*`"Lb"`)
- `"R"`: cable resistance (default: `"len"`*`"Rb"`)


```julia
env.nc.parameters["cable"]
```


    3-element Vector{Any}:
     Dict{Any, Any}("Cb" => 4.0e-7, "Lb" => 0.000264, "Rb" => 0.722, "C" => 0.0004, "i_limit" => 1.0e13, "len" => 1.0, "L" => 0.00025, "R" => 0.208)
     Dict{Any, Any}("Cb" => 4.0e-7, "Lb" => 0.000264, "Rb" => 0.722, "C" => 0.0004, "i_limit" => 1.0e13, "len" => 1.0, "L" => 0.00025, "R" => 0.208)
     Dict{Any, Any}("Cb" => 4.0e-7, "Lb" => 0.000264, "Rb" => 0.722, "C" => 0.0004, "i_limit" => 1.0e13, "len" => 1.0, "L" => 0.00025, "R" => 0.208)


For all information about the parameters which can be defined investigate the parameter dict itself or have a look into the documentation (`LINK`).

The next step would be to interact with the env to run an experiment.
To learn how to do this, see the `Env_Interaction_DEMO.ipynb`.
