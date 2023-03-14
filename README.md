# RANK Example for the BASE for HANK toolbox.

This package provides the example of a RANK economy for the [BASEforHANK](https://github.com/BASEforHANK/BASEtoolbox.jl/) Toolbox.

The module runs with Julia 1.8.2. We recommend to use [Julia for VSCode IDE](https://www.julia-vscode.org) as a front-end to Julia. To get started with the example, simply download or clone the folder, e.g. via `git clone`, `cd` to the project directory and call
```julia-repl
(v1.8) pkg> activate .

(BASEtoolbox_RANK) pkg> instantiate
```
This will install all needed packages. For more on Julia environments, see [`Pkg.jl`](https://julialang.github.io/Pkg.jl/v1/environments/#Using-someone-else's-project).

!!! warning
    Before you activate the environment, make sure that you are in the main directory, in which the `Manifest.toml` and `Project.toml` files are located. In case you accidentally activated the environment in a subfolder, empty `.toml` files will be created that you need to delete before proceeding in the correct folder.


The full documentation of the BASEforHANK toolbox can be found [here](https://baseforhank.github.io/BASEtoolbox.jl/).
