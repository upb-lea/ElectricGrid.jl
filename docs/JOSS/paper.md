---
title: 'ElectricGrid.jl - A Julia-based modeling and simulation tool for power electronics-driven electric energy grids'
tags:
  - Julia
  - Electric Grids
  - Microgrids
  - Reinforcement Learning
  - Energy Systems
  - Simulation
  - Testing
  - Control
authors:
  - name: Oliver Wallscheid
    orcid: 0000-0001-9362-8777
    affiliation: 1
  - name: Sebastian Peitz
    orcid: 0000-0002-3389-793X
    affiliation: 2
  - name: Jan Stenner
    affiliation: 2
  - name: Daniel Weber
    orcid: 0000-0003-3367-5998 
    affiliation: 1
  - name: Septimus Boshoff
    affiliation: 1
  - name: Marvin Meyer
    affiliation: 1
  - name: Vikas Chidananda
    affiliation: 2
  - name: Oliver Schweins
    affiliation: 1
affiliations:
  - name: Chair of Power Electronics and Electrical Drives, Paderborn University, Paderborn, Germany
    index: 1
  - name: Chair of Data Science for Engineering, Paderborn University, Paderborn, Germany
    index: 2
date: 28 April 2023
bibliography: paper.bib
---

# Summary

The ElectricGrid.jl toolbox provides a transient simulation framework
for electric energy grids based on power electronic converters. 
With a few lines of code, a parameterised electric grid model can 
be initialised in Julia. 
An example grid is shown in the figure below.

![Exemplary electric energy grid in a simplified single phase representation.\label{fig:grid}](ExampleGrid.png)

By means of  state-space models (set of first order ordinary differential equations), electrical energy grids can be constructed and simulated in a short time, which can be utilised
for synthetic data generation (e.g., for data-driven meta-modelling) and to compare grid control methods.


# Statement of need

Decentralised, electrical energy networks have special demands on operating and control procedures 
to ensure a continuous and efficient energy supply and
simultaneously play an important role in the integration of renewable energy sources [@Guerrero2013]. This applies both in connection with conventional centralised power grids and for islanded microgrids in remote areas [@Lund2017].
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
which cannot be evaluated a priori [@Garcia2015].
Besides RL methods, being data-driven, result in a model-free and self-adaptive controller design with little human effort labeling them a promising tool for controlling unknown or changing systems targeting the above described challanges.
However, there is a need for further research into the requirements for energy networks in terms of safety, robustness and availability before RL-based controllers can be used in real applications [@ZhangRL2018], [@GlavicRL2018].
 

``ElectricGrid.jl`` is a Julia package for setting up realistic electric grid simulations with support for control options. A number of parameters are made avaible to the user to evaluate the various control options. If no details are given, all parameters are generated automatically, either through randomness on a physically meaningful basis or by verified design methods.
This enables both experts from the field of electrical energy networks to test certain configurations and 
experts from the field of artificial intelligence to test new control approaches
without any prior knowledge of electrical engineering.
Therefore, ``ElectricGrid.jl`` is designed to be used by students, academics, and industrial researchers in the field of simulation and data-driven analysis of electrical energy systems. 
The primary objective of the toolbox is to facilitate entry for new users into 
the modeling, control, and testing of
small to large scale electric power grids and to provide a platform on which different control methods (including RL) can be compared under defined conditions (benchmarks).

The experiments are based on dynamic simulations in the time domain which allows for accurate control
and test investigations during transients and steady state. 
This is an essential difference to already available open-
source solutions for the simulation of electrical energy grids like PyPSA [@Brown2017], Powermodels.jl [@Coffrin2018] and pandapower 
[@Thurner2018] which, in contrast, usually perform the calculations 
in a (quasi)-stationary state. Also these frameworks tend to focus on large-scale power systems at the transmission and distribution grid level, 
which does not allow an evaluation of, e.g., control on component level in case of load fluctuations in the grid.
In addition, few tools like [@lara2023revisiting] already exist in Julia, which offer dynamic simulations. 
However, the latter is based on different simplifications (e.g. assumption of a symmetric grid, fixed frequency, ...) and 
also do not offer an interface to RL toolboxes. In contrast to this, to ensure seamless integration of the control algorithms, 
the Gymnasium-based API [@FaramaFoundation2023] is implemented in ``ElectricGrid.jl``. For other projects implementing the Gymnasium-based API, such as GridAlive (based on Grid2Op [@Donnot2020]) and ChroniX2Grid [@Marot2020], the focus remains on top-down control and steady-state models.
Therefore, this contribution provides a tool to close this gap. 


# Interfaces for control and reinforcement learning  


The API is designed to provide a user-friendly interface to connect a modeled electric energy grid 
with a wide range of classical control methods like shown in the figure below.

![Overview of the functionality and interconnections of the ElectricGrid.jl framework.\label{fig:jeg}](Overview_EG.png)
 
Already provided are classic controllers (i.e., industry standard contollers) like linear feedback proportional integral (PI) in the direct-quadrature-zero (DQ0) rotating reference frame. These control methods can be used out of the box including automatical tuning procedures. 
Many basic auxiliary functionalities for the essential operation of electric power grids are provided too such
as coordinate transformations for basic controller classes, data logging, measurement of real and imaginary powers, 
and phase-locked loops for frequency and phase angle extraction. 
The interface provided by @Tian2020Reinforcement is also available for training 
data-driven control approaches like RL.
This enables users who want to integrate contemporary open-source Julia-based RL toolboxes
such as ``ReinforcementLearning.jl`` [@Tian2020Reinforcement].
Following this structure, nearly every control approach, including data-driven RL, can 
be implemented and tested with ``ElectricGrid.jl`` in a relatively short amount of time. 



# Features

The ``ElectricGrid.jl`` toolbox provides the following key features:


* Framework to set up an experiment with a parameterised energy grid in a few lines of code. 

* Dynamic simulation of electricity grids on component level including single and multi-phase systems as well as AC and DC operation. 

* Calculation, evaluation and logging of every single time step covering states, action and auxiliary quantities. 

* Large variety of predefined and parameterisable controllers (droop, VSG, swing, active-reactive) are available.

* Interesting use cases applying data-driven learning.

# Examples
For illustration and interactive introduction, Jupyter Notebooks are available for each topic.
These provide clear and easy-to-expand examples of: 
 - [Utilising ElectricGrid.jl to build an energy grid](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/Env_Create_DEMO.ipynb),
 - [Theoretical principles behind the calculations](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/NodeConstructor_Theory_DEMO.ipynb),
 - [Applying classic controllers on the electrical grid](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/3_Classical_Controllers_Droop.ipynb),
 - [Training an RL agent on the electrical grid](https://github.com/upb-lea/ElectricGrid.jl/blob/main/examples/notebooks/RL_Single_Agent_DEMO.ipynb).


# Availability and installation

``ElectricGrid.jl`` is supported and tested on Linux, Windows and macOS.
The package should be installed using the Julia package manager. In a Julia terminal run the follwing:
```
]
add ElectricGrid
```

Alternatively it can also be installed from the Github source code. To do that, clone the repository, start Julia, activate the project by pressing `]` to access Pkg mode and then `activate path/to/ElectricGrid` or `activate .` If you started Julia in your ElectricGrid directory and afterwards run `instantiate`.
  
The source code, guide and 
examples are available on the GitHub repository (https://github.com/upb-lea/JuliaElectricGrid.jl). 

# Individual contributions of the authors

Following are shown the main fields of each individual contributor of ElectricGrid.jl: 

* O. Wallscheid: Concept design and idea generation, testing and technical feedback, administrative project management

* S. Peitz: Administrative project management, concept-oriented feedback

* J. Stenner: API RL framework, API environment framework, basic system architecture

* D. Weber: Application examples, API environment framework, basic system architecture, unit tests

* S. Boshoff: Application examples, primary controllers in DQ0 frame, decentralised secondary controllers, Luenberger observers, stochastic processes, inverter filter design, cable design, unit tests 

* M. Meyer: Basic system architecture, application examples, unit tests

* V. Chidananda: System analytics, unit tests

* O. Schweins: Basic system architecture




# Acknowledgements

The authors would also like to acknowledge 
the funding and support of this work by the Federal Ministry of Education and Research of Germany (grant number 011S21064).  

# References

