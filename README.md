
# ElectricGrid.jl

<img align="right" width="150" height="200" src="docs/logo.png">

| [**Reference docs**](https://upb-lea.github.io/ElectricGrid.jl/dev/)
| [**Install guide**](#installation)
| [**Quickstart**](#getting-started)
| [**Release notes**](https://github.com/upb-lea/ElectricGrid.jl/releases/new)

[![DOI](https://joss.theoj.org/papers/10.21105/joss.05616/status.svg)](https://doi.org/10.21105/joss.05616)
[![Build Status](https://github.com/upb-lea/ElectricGrid.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/upb-lea/ElectricGrid.jl/actions/workflows/CI.yml)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/upb-lea/ElectricGrid.jl/blob/main/LICENSE)




ElectricGrid.jl is a library for setting up realistic electric grid simulations with extensive support for control options. With ElectricGrid.jl you can
- create a simulation environment for an electric grid by defining its sources, loads, and cable connections,
- set detailed parameters of your electric components or let them be auto-generated,
- choose different control modes for each power electronic converter in your system and
- train RL agents as controllers with the built-in DDPG/TD3 implementations (based on [Flux.jl](https://fluxml.ai/)) or write your own. The environment and agent interface follows the design of [ReinforcementLearning.jl](https://juliareinforcementlearning.org/) v0.10, which ElectricGrid now bundles internally instead of depending on it.


![ElectricGrid Framework](docs/src/assets/Overview_EG.png)

## Installation

ElectricGrid.jl is tested with Julia 1.11 and 1.12.

- Installation using the Julia package manager (recommended if you want to use ElectricGrid in your project):
  - In a Julia terminal run the following:
```
import Pkg
Pkg.add("ElectricGrid")
```
or press `]` in the Julia Repl to enter Pkg mode and then run
```
add ElectricGrid
```

- Install from GitHub source (recommended if you want to run the example notebooks and scripts):
  - Clone the git and navigate to the directory
```
git clone https://github.com/upb-lea/ElectricGrid.jl.git
```

## Getting Started

To get started with ElectricGrid.jl the following interactive notebooks are useful. They show how to use the ElectricGrid.jl framework to build and simulate the dynamics of an electric power grid controlled via classic controllers or train common RL agents for different control tasks:
* [Create an environment with ElectricGrid.jl](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/Env_Create.ipynb)
* [Theory behind ElectricGrid.jl - Modelling Dynamics using Linear State-Space Systems](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/NodeConstructor_Theory.ipynb)
* [Classic Controlled Electric Power Grids - State-of-the-Art](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/Classical_Controllers_1_Swing.ipynb)
* [Use RL Agents in the ElectricGrid.jl Framework](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/RL_Single_Agent.ipynb)

An overview of all parameters defining the experiment setting in regard to the electric grid can be found here:
* [Default Parameters](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/Default_Parameters.ipynb)


To run a simple example, the following few lines of code can be executed:

```
using ElectricGrid

env =  ElectricGridEnv(num_sources = 1, num_loads = 1)
Multi_Agent =  SetupAgents(env)
hook =  Simulate(Multi_Agent, env)
RenderHookResults(hook = hook)
```

This is a minimal example of a full ElectricGrid.jl setup. 
There should also appear a plot that looks like this:
![output of the minimal example](docs/src/assets/output1.png)


## Using the GUI

The current version of ElectricGrid features a graphical user interface (GUI) that helps with setting up a simulation.
This is built on the library [QML.jl](https://github.com/JuliaGraphics/QML.jl), which installs normally from the Julia package registry (it is a dependency of ElectricGrid) and works with the current Qt6-based QML.jl v0.13.
To start it, run `julia --project=. gui/ElectricGridGUI.jl` from the repository root.

![GUI example](docs/src/assets/gui_example.png)

Usage of the GUI is explained in the [GUI section in the docs](https://upb-lea.github.io/ElectricGrid.jl/dev/Gui/).
