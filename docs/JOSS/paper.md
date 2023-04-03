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
  - name: Septimus Boshoff
    affiliation: 1
  - name: Marvin Meyer
    affiliation: 1
  - Vikas Chidananda
    affiliation: 2
  - name: Oliver Schweins
    affiliation: 1
  
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
With a few lines of code, a parameterized energy grid can 
be initialized in Julia. 
An example is shown in the figure below.

![\label{fig:grid}](ExampleGrid.png)
_Fig. 1:  Exemplary electric energy grid in single phase representation._

By means of a linear time invariant state-space system, the system is evolved, which can then be
for synthetic data generation to test novel control methods and compare them 
with state-of-the-art controllers.


# Background on electric grids and their control


Decentralized, electrical energy networks have special demands on operating and control procedures 
to ensure a continuous and efficient energy supply and
simultaneously play an important role in the integration of renewable energy sources. [@Guerrero2013] 
This applies both in connection with conventional power grids and for power supply in remote areas [@Lund2017].
Due to their high efficiency and flexibility, the integration of power electronic converters 
in electrical energy grids is increasing. 
Power electronics describes the application of solid-state electronics to the control and
conversion of electric power, which is largely performed with semiconductor switching 
devices such as diodes or power transistors.
This includes energy conversion in terms of voltage and current amplitude, frequency
and phase angle, as well as the number of phases between two or more
electrical energy systems to be connected.


Controlling (decentralized) electric grids is a challenging task due to their stochastic, heterogeneous 
and volatile characteristics.
At the same time, high requirements are made with regard to aspects such as safety, quality and availability.
This results in a high demand for comprehensive testing of new control concepts during their development phase and comparisons with the state
of the art to ensure their feasibility.
This applies in particular to data-driven control approaches such as 
reinforcement learning (RL), the stability and operating behavior of
which cannot be evaluated a priori [@Garcia2015].


# State of field 

``JEG`` is a Julia package for setting up realistic electric grid simulations with support for control options. A number of parameters are made avaible to the user to evaluate the various control options. If no details are given, all parameters are generated automatically, either through randomness on a physically meaningful basis or by verified design methods.
This enables both experts from the field of electrical energy networks to test certain configurations and 
experts from the field of artificial intelligence to test new control approaches
without any prior knowledge of electrical engineering.
The experiments are based on dynamic simulations in the time domain which allows for accurate control
and test investigations during transients and steady-state. 
This is an essential difference to already available open-
source solutions for the simulation of electrical energy grids, which, in contrast, usually perform the calculations 
in a (quasi)-stationar state [@Coffrin2018],
which does not allow an evaluation of e.g. control on component level in case of load fluctuations in the grid.




# Interfaces for control and reinforcement learning  


The API is designed to provide a user-friendly interface to connect a modeled electric energy grid 
with a wide range of classical control methods like shown in the figure below.

![\label{fig:jeg}](OverviewJEG.png)
_Fig. 2:  Overview of the functionality and interconnections of the JEG framework._
 
Already provided are classic controllers (i.e. industry standard contollers) like linear feedback Proportional Integral (PI) in the Direct-Quadrature-Zero (DQ0) rotating reference frame. These control methods can be used out of the box. That is, the controllers are automatically tuned with Luenberger observers. 
Many basic auxiliary functionalities for the essential operation of electric power grids are provided by JEG such
as coordinate transformations for basic controller classes, data logging, measurement of real and imaginary powers, 
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


* Framework to set up an experiment with a parameterized energy grid in a few lines of code, 

* Dynamic simulation of local electricity grids on component level including single and multi-phase systems as well as AC and DC operation. 

* Calculation, evaluation and logging of every single time step covering states, action and auxiliary quantities. 

* Large variety of predefined and parameterizable controllers (droop, VSG, swing, active-reactive) are available.

* Interesting use cases applying data-driven learning.

# Examples
For illustration and interactive introduction, jupyter notebooks are available for each topic.
These provide clear and easy-to-expand examples of: 
 - [Utilising JEG to build an energy grid](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/Env_Create_DEMO.ipynb)
 - [Theoretical principles behind the calculations](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/NodeConstructor_Theory_DEMO.ipynb)
 - [Applying classic controllers on the electrical grid](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/1_Auxiliaries_OU_process.ipynb) and
 - [Training an RL agent on the electrical grid](https://github.com/upb-lea/JuliaElectricGrid.jl/blob/main/examples/notebooks/NodeConstructor_Theory_DEMO.ipynb).


# Availability and installation

``JEG`` is supported and tested on Linux and Windows. 
The package should be installed ... @Jan -> #43: https://github.com/upb-lea/JuliaElectricGrid.jl/issues/43
The source code, guide and 
examples are available on the GitHub repository (https://github.com/upb-lea/JuliaElectricGrid.jl). 

# Individual contributions of the authors

Following are shown the main fields of each individual contributor of JEG: 

* O. Wallscheid: Concept design and idea generation, testing and technical feedback, administrative project management

* S. Peitz: Administrative project management, concept-oriented feedback

* J. Stenner: API RL framework

* D. Weber: Application examples, API environment framework, unit tests

* S. Boshoff: Application examples, primary controllers in DQ0 frame, decentralised secondary controllers, Luenberger observers, stochastic processes, inverter filter design, cable design, unit tests 

* M. Meyer: Basic system architecture, application examples, unit tests

* V. Chidananda:

* O. Schweins: Basic system architecture




# Acknowledgements

The authors would also like to acknowledge 
the funding and support of this work by the BMBF research grant. 

# References

