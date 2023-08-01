# ElectricGrid.jl

ElectricGrid.jl is a library for setting up realistic electric grid simulations with extensive support for control options. With ElectricGrid.jl you can
- Create a simulation environment for a electric grid by defining its sources, loads and cable connections
- Set detailled parameters of your electric components - or let them be auto-generated
- Choose different control modes for each source in your system
- Use the agent architecture of ReinforcementLearning.jl to either train RL agents as controllers or write your own ones

The getting started pages of this documentation will guide you through the principles of working with ElectricGrid.jl by following different examples and explaining them step by step. The example scripts can also be found in the examples folder of the ElectricGrid.jl repository.

## Installation
- Installation using the julia package manager (recommended if you want to use ElectricGrid in your project):
In a julia terminal run the follwing:
```
import Pkg
Pkg.add("ElectricGrid")
```

- Install from Github source (recommended if you want to run the example notebooks and scripts):
  - Clone the git and navigate to the directory
```
git clone https://github.com/upb-lea/ElectricGrid.jl.git
```





