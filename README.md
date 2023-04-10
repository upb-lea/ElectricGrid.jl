
# Julia Electric Grid - JEG  

<img align="right" width="150" height="200" src="docs/logo.png">

| [**Reference docs**](https://upb-lea.github.io/JuliaElectricGrid.jl/dev/)
| [**Install guide**](#installation)
| [**Quickstart**](#getting-started)
| [**Release notes**](https://github.com/upb-lea/JuliaElectricGrid.jl/releases/new)

[![Build Status](https://github.com/upb-lea/JuliaElectricGrid.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/upb-lea/JuliaElectricGrid.jl/actions/workflows/CI.yml)
[![License](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/LICENSE)




JuliaElectricGrid, or JEG for short, is a library for setting up realistic electric grid simulations with extensive support for control options. With JEG you can,
- create a simulation environment for an electric grid by defining its sources, loads, and cable connections
- set detailed parameters of your electric components - or let them be auto-generated
- choose different control modes for each source in your system
- use the agent architecture of [ReinforcementLearning.jl](https://juliareinforcementlearning.org/) to either train RL agents as controllers or write your own ones


![JEG Framework](docs/src/assets/OverviewJEG.png)

## Installation
- Installation using the julia package manager (recommended):
In a julia terminal run the follwing:
```
]
activate .
add JEG
```

- Install from Github source:
  - Clone the git and navigate to the directory
```
git clone https://github.com/upb-lea/JuliaElectricGrid.jl.git
```
  - activate Julia

  - activate the project by pressing `]`to access pkg mode and then `activate path/to/JEG` or `activate .` if you started julia in your JuliaElectricGrid directory
  - run `instantiate`

## Getting Started

To get started with JEG the following interactive notebooks are useful. They show how to use the JEG framework to build and simulate the dynamics of an electric power grid controlled via classic controllers or train common RL agents for different control tasks:
* [Create an environment with JEG](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/Env_Create_DEMO.ipynb)
* [Theroy behind JEG - Modelling Dynamics using Linear State-Space Systems](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/NodeConstructor_Theory_DEMO.ipynb)
* [Classic Controlled Electric Power Grids - State of the Art](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/Classical_Controllers_Introduction.ipynb)
* [Use RL Agents in the JEG Framework](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/RL_Single_Agent_DEMO.ipynb)

An overview of all parameters defining the experiment setting with regards to the electric grid can be found here:
* [Default Parameters](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/Default_Parameters.ipynb)


To run a simple example, the following few lines of code can be executed:

```
using JEG

env =  ElectricGridEnv(num_sources = 1, num_loads = 1)
Multi_Agent =  SetupAgents(env)
hook =  Simulate(Multi_Agent, env)
RenderHookResults(hook = hook)
```

This is a minimal example of a full JuliaElectricGrid setup. 
There should also appear a plot that looks like this:
![output of the minimal example](docs/src/assets/output1.png)

