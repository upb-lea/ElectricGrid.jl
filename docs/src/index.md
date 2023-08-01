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

## Statement of need

Decentralised, electrical energy networks have special demands on operating and control procedures 
to ensure a continuous and efficient energy supply and
simultaneously play an important role in the integration of renewable energy sources ([@Guerrero2013](https://ieeexplore.ieee.org/document/6184305)). This applies both in connection with conventional centralised power grids and for islanded microgrids in remote areas ([@Lund2017](https://www.sciencedirect.com/science/article/abs/pii/S0360544217308812)).
Due to their high efficiency and flexibility, power electronic converters have become the standard tool for integrating renewable energy sources, energy storages and loads in 
electrical energy grids. 
The field of power electronics covers the application of solid-state electronics to the control and
conversion of electric power, which is performed with semiconductor switching 
devices such as diodes or power transistors.
This includes energy conversion in terms of voltage and current amplitude, frequency
and phase angle, as well as the number of phases between two or more electrical energy systems to be connected.

Controlling (decentralised) electric grids is a challenging task due to their stochastic, heterogeneous 
and volatile characteristics (in particular regarding the connected loads).
At the same time, high requirements are made with regard to aspects such as safety, quality and availability.
This results in a high demand for comprehensive testing of new control concepts during their development phase and comparisons with the state
of the art to ensure their feasibility.
This applies in particular to emerging data-driven control approaches such as 
reinforcement learning (RL), the stability and operating behavior of
which cannot be evaluated a priori ([@Garcia2015](https://www.jmlr.org/papers/volume16/garcia15a/garcia15a.pdf)).
Besides RL methods, being data-driven, result in a model-free and self-adaptive controller design with little human effort labeling them a promising tool for controlling unknown or changing systems targeting the above described challanges.
However, there is a need for further research into the requirements for energy networks in terms of safety, robustness and availability before RL-based controllers can be used in real applications ([@ZhangRL2018](https://ieeexplore.ieee.org/abstract/document/8468674), [@GlavicRL2018](https://www.sciencedirect.com/science/article/pii/S2405896317317238)).


``ElectricGrid.jl`` is a Julia package for setting up realistic electric grid simulations with support for control options. A number of parameters are made avaible to the user to evaluate the various control options. If no details are given, all parameters are generated automatically, either through randomness on a physically meaningful basis or by verified design methods.
This enables both experts from the field of electrical energy networks to test certain configurations and 
experts from the field of artificial intelligence to test new control approaches
without any prior knowledge of electrical engineering.
Therefore, ``ElectricGrid.jl`` is designed to be used by students, academics, and industrial researchers in the field of simulation and data-driven analysis of electrical energy systems. 
The primary objective of the toolbox is to facilitate entry for new users into 
the modeling, control, and testing of
small to large scale electric power grids and to provide a platform on which different control methods (including RL) can be compared under defined conditions (benchmarks).

The experiments are based on dynamic simulations in the time domain which allows for accurate control
and test investigations during transients and steady state down to component level. 
This is an essential difference to already available open-source solutions for the simulation of electrical energy grids like PyPSA ([@Brown2017](https://pypsa.readthedocs.io/en/latest/)), Powermodels.jl ([@Coffrin2018](https://arxiv.org/abs/1711.01728)) and pandapower 
([@Thurner2018](https://arxiv.org/pdf/1709.06743.pdf)) which, in contrast, usually perform the calculations 
in a (quasi)-stationary state. Also these frameworks tend to focus on large-scale power systems at the transmission and distribution grid level, 
which does not allow an evaluation of, e.g., control on component level in case of load fluctuations in the grid.
In addition, few tools like [@lara2023revisiting](https://arxiv.org/abs/2301.10043) already exist in Julia, which offer dynamic simulations. 
However, the latter is based on different simplifications (e.g., assumption of a symmetric grid, fixed frequency, ...) and 
also do not offer an interface to RL toolboxes.
To ensure a seamless integration of the control algorithms, the Gymnasium-based API ([@FaramaFoundation2023](https://github.com/Farama-Foundation/Gymnasium)) should be used, which has been established as a standard in recent years.
For other projects implementing the Gymnasium-based API, such as GridAlive (based on Grid2Op [@Donnot2020](https://github.com/rte-france/Grid2Op)) and ChroniX2Grid ([@Marot2020](https://github.com/BDonnot/ChroniX2Grid)), the focus remains on top-down control and steady-state models.
Therefore, ``ElectricGrid.jl`` provides a tool to close these highlighted gaps.


