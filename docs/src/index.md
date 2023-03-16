# JuliaElectricGrid.jl

JuliaElectricGrid, or JEG for short, is a library for setting up realistic electric grid simulations with extensive support for control options. With JEG you can
- Create a simulation environment for a electric grid by defining its sources, loads and cable connections
- Set detailled parameters of your electric components - or let them be auto-generated
- Choose different control modes for each source in your system
- Use the agent architecture of ReinforcementLearning.jl to either train RL agents as controllers or write your own ones

The getting started pages of this documentation will guide you through the principles of working with JEG by following different examples and explaining them step by step. The example scripts can also be found in the examples folder of the JEG repository.

## Installation

For now you have to download or clone the github repository at https://github.com/upb-lea/JuliaElectricGrid.jl and set it up:
- start up julia
- activate the project by pressing `]`to access pkg mode and then `activate path/to/JEG` or `activate .` if you started julia in your JuliaElectricGrid directory
- run `instantiate`

Now you should be ready to run the examples and work with JEG.





