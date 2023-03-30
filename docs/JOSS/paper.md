---
title: 'JEG: Julia Electric Grid'
tags:
  - Julia
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
    orcid: 0009
    affiliation: 2
  - name: Jan Stenner
    affiliation: 2
  - name: Daniel Weber
    orcid: 0000-0003-3367-5998 
    affiliation: 1
  - name: Septimus Boschhoff
    affiliation: 1
  - name: Marvin Meyer
    affiliation: 1
  - SHK/WHB?
  
affiliations:
 - name: Chair of Power Electronics and Electrical Drives, Paderborn University, Paderborn, Germany
   index: 1
- name: Chair of Data Science for Engineering, Paderborn University, Paderborn, Germany
   index: 1
date: XX April 2023
bibliography: paper.bib
---

# Summary

The Julia Electric Grid (JEG) toolbox provides a transient simulation framework
for energy grids based on power electronic converters. 
With a few lines of code, an arbitrarily detailed parameterized energy grid can 
be initialized in Julia. 
An example is shown in the figure below.

![\label{fig:grid}](ExampleGrid.png)
_Fig. 1:  Exemplary electric energy grid in single phase representation._

By means of a state-space system, the system response is determined in the time domain which can 
be used for synthetic data generation to test novel control methods and compare them 
with state-of-the-art controllers.


# Background on electric grids and their control


Decentralized, electrical energy networks place special demands on operating and control procedures to ensure a continuous and efficient energy supply and
at the same time play an important role in the integration of renewable energy sources. 
This applies both in connection with conventional power grids and for power supply in remote areas [@Lund2017]. 
Due to their high efficiency and flexibility, power electronic converters are largely used to
drive modern MSG.
Power electronics describes the application of solid-state electronics to the control and
conversion of electric power, which is largely performed with semiconductor switching 
devices such as diodes or power transistors.
This includes energy conversion in terms of voltage and current amplitude, frequency
and phase angle, as well as the number of phases between two or more
electrical energy systems to be connected.


Controlling (decentralized) electric grids is a challenging task due to their stochastic, heterogeneous 
and volatile characteristics
At the same time, high requirements are made with regard to aspects such as safety, quality and availability.
This results in high demand for comprehensive testing of new control 
concepts during their development phase and comparisons with the state
of the art to ensure their feasibility.
This applies in particular to data-driven control approaches such as 
reinforcement learning (RL), the stability and operating behavior of
which cannot be evaluated a priori [@Garcia2015].


# State of field 

JEG is a Julia package for setting up realistic electric grid simulations with 
extensive support for control options.
The parametrisation can be arbitrarily detailed. 
If no details are given, all parameters are calculated randomly on a physically meaningful basis.
This enables both experts from the field of electrical energy networks to test certain configurations and 
experts from the field of artificial intelligence to test new control approaches
without any prior knowledge of electrical engineering.
The experiments are based on dynamic simulations in the time domain which allows for accurate control
and test investigations during transients and steady-state. 
This is an essential difference to already available open-
source solutions for the simulation of electrical energy grid, which, in contrast, usually perform the calculations 
in a (quasi)-stationar state [@Coffrin2018], 
which does not allow an evaluation of e.g. control on component level in case of load fluctuations in the grid.




# Interfaces for control and reinforcement learning  

RL.jl, classic controllers, PLL, transforms, mulicontroller

The API is designed to provide a user-friendly interface to connect a modeled electric energy grid 
with a wide range of classic control methods like shown in the figure below.

![\label{fig:jeg}](OverviewJEG.png)
_Fig. 2:  Overview of the functionality and interconnections of the JEG framework._
 
Already provided classic controllers like linear feedback control methods can be used out of the box or 
be extended by model predictive control techniques. 
Many basic auxiliary functionalities for the essential operation of electric power grids are provided by JEG such
as coordinate transformations for basic controller classes, data logging, 
and phase-locked loops for frequency and phase angle extraction. 
The interface provided by [@Tian2020Reinforcement] is also available for training 
data-driven control approaches like RL.
This enables users who want to integrate contemporary open-source Julia-based RL toolboxes
such as ``ReinforcementLearning.jl`` [@Tian2020Reinforcement].
Following this structure, nearly every control approach, including data-driven RL, can 
be implemented and tested with ``JEG`` in a relatively short amount of time. 


# Intended use and targeted audience

``JEG`` is designed to be used by students, academics, and industrial researchers 
in the field of control and energy engineering and data science. 
The primary objective of the toolbox is to facilitate entry for new users into 
the modeling, control, and testing of
small to large scale electric power grids and to provide a platform on which different control methods (including RL) 
can be compared under defined conditions (benchmarks).



# Features

The ``JEG`` toolbox provides the following key features:


* A library for the scalable and flexible design of local electricity grids in OpenModelica.
Users can select between a wide range of different grid components and connect them in a plug-and-play approach.

* Dynamic simulation of local electricity grids on component level including single and multi-phase systems as well as AC and DC operation. 

* Calculation, evaluation and logging of every single time step covering states, action and auxiliary quantities. 

* Large variety of predefined and parameterizable controllers (droop, VSG, voltage, current) available, easy implementation of user-defined control structures possible.

* Interesting use cases applying data-driven learning.

# Examples
For illustration and interactive introduction, jupyter notebooks are available for each topic.
These provide clear and easy-to-expand examples of: 
 - [how JEG can be used to set up a network](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/Env_Create_DEMO.ipynb),
 - [the theoretical principles behind the calculations](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/NodeConstructor_Theory_DEMO.ipynb),
 - [how the classic controllers work](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/1_Auxiliaries_OU_process.ipynb) and
 - [how to train an RL agent](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/NodeConstructor_Theory_DEMO.ipynb).


# Availability and installation

``JEG`` is supported and tested on Linux and Windows. 
The package should be installed ... @Jan -> #43: https://github.com/upb-lea/JuliaElectricGrid.jl/issues/43
The source code, guide and 
examples are available on the GitHub repository (https://github.com/upb-lea/JuliaElectricGrid.jl). 

# Individual contributions of the authors

Following are shown the main fields of each individual contributor of OMG: 

* O. Wallscheid: Concept design and idea generation, testing and technical feedback, administrative project management

* S. Peitz: Administrative project management, concept-oriented feedback

* J. Stenner: API RL framework

* D. Weber: Application examples, API environment framework, unit tests

* S. Boshoff: Application examples, classic controllers

* M. Meyer: Application examples, theoretical ...

* SHK/WHB?: 




# Acknowledgements

The authors would also like to acknowledge 
the funding and support of this work by the BMBF research grant. 

# References

