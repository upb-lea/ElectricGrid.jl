
![Dare Logo](public/logo.png)

## Instantiating and managing the "dare" julia environment

### Setting up the project with VSCode
It's recommended to work in this project with Visual Studio Code.
- In VSCode, use "Open Folder" and select the "dare" folder that contains the `Project.toml`. This way, VSCode should automatically recognize the julia environment called `dare` - or it should ask you for selecting an environment and `path/to/dare` will be in the list to choose.
- Open a new terminal - it should start in the "dare" folder. Then start julia in it (the command is just `julia` if you included it in your systems PATH variable). Now, type `]` to access the package manager (the prompt color switches from green to blue) and then run `activate .` (with the dot). The project gets activated and it shows the name in the brackets. \
 ![activated project](public/doc/activate.png) \
 Afterwards, run `instantiate`. A progress starts that automatically installs all the needed packages for the julia environment called `dare`.
 ![instantiate project](public/doc/instantiate.png)
 It gets the information about the required libraries from the `Project.toml` and creates a new file called `Manifest.toml` which stores information about the all the sub-dependencies of these libraries installed. `Manifest.toml` might be different from machine to machine and is therefore a local file that is also excluded from the git via the `.gitignore` file. 
- Now, if you "Execute active File in REPL" it should automatically activate the project before running the file.
 ![run in REPL](public/doc/RunInREPL.png)
 It might take some time to precompile packages (which is only done once as long as package versions don't change) and also to precompile parts of your code (which has to be done any time you run in a new REPL).