
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
- create a simulation environment for an electric grid by defining its sources, loads, and cable connections
- set detailed parameters of your electric components - or let them be auto-generated
- choose different control modes for each source in your system
- use the agent architecture of [ReinforcementLearning.jl](https://juliareinforcementlearning.org/) to either train RL agents as controllers or write your own ones


![ElectricGrid Framework](docs/src/assets/Overview_EG.png)

## Installation
- Installation using the julia package manager (recommended if you want to use ElectricGrid in your project):
  - In a Julia terminal run the follwing:
```
import Pkg
Pkg.add("ElectricGrid")
```
or press `]` in the Julia Repl to enter Pkg mode and then run
```
add ElectricGrid
```

- Install from Github source (recommended if you want to run the example notebooks and scripts):
  - Clone the git and navigate to the directory
```
git clone https://github.com/upb-lea/ElectricGrid.jl.git
```

## Getting Started

To get started with ElectricGrid.jl the following interactive notebooks are useful. They show how to use the ElectricGrid.jl framework to build and simulate the dynamics of an electric power grid controlled via classic controllers or train common RL agents for different control tasks:
* [Create an environment with ElectricGrid.jl](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/Env_Create_DEMO.ipynb)
* [Theory behind ElectricGrid.jl - Modelling Dynamics using Linear State-Space Systems](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/NodeConstructor_Theory_DEMO.ipynb)
* [Classic Controlled Electric Power Grids - State of the Art](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/Classical_Controllers_Introduction.ipynb)
* [Use RL Agents in the ElectricGrid.jl Framework](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/RL_Single_Agent_DEMO.ipynb)

An overview of all parameters defining the experiment setting with regards to the electric grid can be found here:
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
This is built on the library [QML.jl](https://github.com/JuliaGraphics/QML.jl), that, at the time of writing, stopped working in it's current release version.
For that reason it is **required to clone this codebase** and install `QML.jl` in it's github main state manually if you want to use the gui.

```
import Pkg
Pkg.add("QML#main")
```
or press `]` in the Julia Repl to enter Pkg mode and then run
```
add QML#main
```

![gui example](docs/src/assets/gui_example.png)

Usage of the GUI is explained in the [GUI section in the docs](https://upb-lea.github.io/ElectricGrid.jl/dev/Gui/).
