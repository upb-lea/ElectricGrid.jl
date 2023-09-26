# Graphical User Interface

![gui example](docs/src/assets/gui_example.png)


## Preliminary steps for using the GUI

The current version of ElectricGrid features a graphical user interface (GUI) that helps with setting up a simulation.
This is built on the library [QML.jl](https://github.com/JuliaGraphics/QML.jl), that, at the time of writing, stopped working in it's current release version.
For that reason it is required to install `QML.jl` in it's github main state manually if you want to use the gui.

```
import Pkg
Pkg.add("QML#main")
```
or press `]` in the Julia Repl to enter Pkg mode and then run
```
add QML#main
```


## Usage

To start the GUI simply run the file `gui\ElectricGridGUI.jl`. If everything works well a blank white window will open.

### Adding Sources and/or Loads

