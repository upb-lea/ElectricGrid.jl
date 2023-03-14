# NodeConstructor Demo

The NodeConstructor has the goal to convert the physical system, here the grid, into a mathematical model. This mathematical model can then be used to simulate the behavior of the network. This notebook addresses the following points:

- ### Introduction to the NodeConstructor
- ### Representation of the physical grid
- ### Bulding of the state-space matrices

![Illustration of where the NodeConstructor is located](./assets/OverviewJEG.png)

The NodeConstructor is a part of the enviroment (see graphic). Its function is to create the mathematical model of the grid to be simulated. The inputs for the NodeConstructor are the specifications of the grid, which will be discussed in more detail subsequently. Within the NodeConstructor, an ODE system is then created based on the underlying physical models and properties of the individual components. Based on this ODE system, the necessary matrices can then be extracted.

In the following, we will discuss how a grid can be transferred to the NodeConstructor.







```@autodocs
Modules = [JEG]
Order   = [:function, :type]
```