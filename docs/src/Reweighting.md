# [Simulation Reweighting](@id simulation_reweighting)

!!! warning
    This is an experimental feature. Breaking changes may occur without 
    a breaking package release.

Computes the new weight for a frame of a simulation based on the energy difference between the perturbed and non-perturbed original sampling

This resource is based on the Free Energy Perturbation Theory (FEP) in the Molecular Dynamics context. Most of the time, each frame will contribute equally
for calculations of some thermodynamic property, however, we can apply a perturbation on one or multiple types of atomic
interactions (e.g. making water oxygen and protein carbonyl oxygen interaction more repulsive), making these frames to have different normalized statistical contributions so that we can 
possibly preview the outcome of a new simulation with these modifications.

## How to use it
```julia-repl
julia> using MolSimToolkit
```

## Setting initial parameters
Firstly, we define the ```simulation``` object and set the atoms that will determine which interactions will be perturbed:

```julia-repl
julia> using MolSimToolkit, PDBTools

julia> testdir = "$(@__DIR__)/test"
"/home/lucasv/.julia/dev/MolSimToolkit/src/Reweighting/test"

julia> simulation = Simulation("$testdir/Testing_reweighting.pdb", "/$testdir/Testing_reweighting_10_frames_trajectory.xtc")
Simulation 
    Atom type: PDBTools.Atom
    PDB file: /home/lucasv/.julia/dev/MolSimToolkit/src/Resampling/test/Testing_resampling.pdb
    Trajectory file: /home/lucasv/.julia/dev/MolSimToolkit/src/Resampling/test/Testing_reweighting_10_frames_trajectory.xtc
    Total number of frames: 10
    Frame range: 1:1:10
    Number of frames in range: 10
    Current frame: nothing

julia> i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

julia> i2 = PDBTools.selindex(atoms(simulation), "protein and name O")
```
In the picture below, it is shown the interaction which will be perturbed in this system:

!!!!INSERT FIGURE

## Setting perturbation function
In order to obtain these weights, we have to use two functions: the ```reweight``` function, which will calculate each weight and the ```perturbation``` function, responsible for taking each computed distance between atomic pairs in every frame and determine the resulting energy using these distances in that particular frame based on the applied perturbation.

So, secondly, we define some "perturbation" function (here we call it ```polynomial_decay```) and set up its parameters. Please, take a look at the interface:

```julia-repl
julia> poly_decay_perturbation(r, α, cut) = α * ((r/cut)^2 - 1)^2
gaussian_decay (generic function with 1 method)

julia> α = 5.e-3
0.005

julia> cut = 10.0
10.0
```

As it can be seen, the function has to receive three parameters: `r` which corresponds to the distance between two selected atoms and some parameter to account a modification and change its magnitude, here, we inserted two of them in the same function, `α`, to change the maximum value of the curve (at r = 1) and `cut`, the distance `r` where the function equals zero. In the image below, we can see the curve and how it changes with different values of `α` and `cut`:

!!! INSERT FIGURE

## Computing the new weights
Finally, using the ```reweight``` function, we pass both the ```simulation``` and the last function anonymously in the input. Again, watch the interface:

```julia-repl

julia> weights = reweight(simulation, r -> gaussian_decay(r, α, cut), i1, i2; cutoff = cut)
```

`r`: the distance between the twos atoms 

`cutoff`: the maximum distance that will be computed between two atoms. The default value is `12.0 Å`.

!!! warning
    It is highly recommended to set the same value of `cutoff` for both `perturbation` and `reweight` functions.
    With this in mind, calculations will be done more quickly and you do not need to worry about your input function
    behaviour above the `cutoff` value, since distances out of the perturbation range will not be computed.

Once the calculations are finished, the resulted interface is shown, like the example below:
```julia-repl
-------------
FRAME WEIGHTS
-------------

Average probability = 0.09999999999999999
standard deviation = 0.01734935311311546

-------------------------------------------
FRAME WEIGHTS RELATIVE TO THE ORIGINAL ONES
-------------------------------------------

Average probability = 0.45655722352062866
standard deviation = 0.0792097248720297

----------------------------------
COMPUTED ENERGY AFTER PERTURBATION
----------------------------------

Average energy = 0.7973733879299723
standard deviation = 0.17177116838361012
```

The data in ```weights``` structure is organized as it follows:

```julia
struct ReweightResults{T<:Real}
    probability::Vector{T}
    relative_probability::Vector{T}
    energy::Vector{T}
end
```

As an example, if we want the absolute weights computed for our simulation:

```julia-repl
julia> weights.probability
10-element Vector{Float64}:
 0.08987791339898044
 0.07326337222373071
 0.0973116226496827
 0.10965810145525891
 0.09829891590498603
 0.0916792371461855
 0.08548699059703141
 0.12480704633057726
 0.09973413264337352
 0.12988266765019355
```

!!! tip
    Note that these values (and, consequently, calculations that use them) are functions of `r`, 
    so in other to avoid mathematical complications, a good piece of advice is to create functions that 
    are continous in the closed interval from zero to the cutoff value.  

## Reference Functions
```@autodocs
Modules = [MolSimToolkit.Reweighting]
Order = [:function, :type]
```