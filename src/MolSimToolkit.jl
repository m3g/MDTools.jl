module MolSimToolkit

import Chemfiles
import PDBTools
import OffsetArrays
import LaTeXStrings # only because Aqua complains: used in the Plotting extensions

using TestItems: @testitem
using StaticArrays: FieldVector, SMatrix, MVector
using LinearAlgebra: norm, cross, dot
using Reexport: @reexport
using ProgressMeter: Progress, next!, @showprogress
using Statistics: mean

export wrap, wrap_to_first
export distances
export dihedral, dihedrals, average_dihedrals
export align, align!, rmsd, rmsd_matrix
export intermittent_correlation
export bulk_coordination
export coordination_number

# Reexported from ProteinSecondaryStructures for convenience
using ProteinSecondaryStructures: dssp_run, stride_run, 
    ss_code, ss_number, ss_name
export dssp_run, stride_run, ss_code, ss_number, ss_name
# SS trajectory functions
export ss_map, ss_mean, ss_heatmap

# New method added to this function, which is reexported
import PDBTools.center_of_mass
export center_of_mass

# Version of the package: used for printing in some places
const version = pkgversion(@__MODULE__)

# Minimal AtomType interface
atomic_mass(atom::PDBTools.Atom) = PDBTools.mass(atom)

# Testing module
include("../test/Testing.jl")

# Data structures
include("./datastructures/Simulation.jl")
include("./datastructures/Positions.jl")

# Coordinate PBC wrapping functions
include("./wrap.jl")

# Miscelaneous functions
include("./miscelaneous_functions/distances.jl")
include("./miscelaneous_functions/dihedrals.jl")
include("./miscelaneous_functions/center_of_mass.jl")
include("./miscelaneous_functions/intermittent_correlation.jl")
include("./miscelaneous_functions/bulk_coordination.jl")
include("./miscelaneous_functions/most_representative_structure.jl")
include("./miscelaneous_functions/coordination_number.jl")

#  Structural alignment
include("./procrustes.jl")

# Analysis functions and modules
include("./BlockAverages.jl")
@reexport using .BlockAverages
include("./MolecularMinimumDistances/MolecularMinimumDistances.jl")
@reexport using .MolecularMinimumDistances
include("./gromacs/remd.jl")
include("./Reweighting/Reweighting.jl")
@reexport using .Reweighting
include("./secondary_structure/secondary_structure.jl")

# Simulation setup facilities
include("./PackmolInputCreator/PackmolInputCreator.jl")

# Cluster managemeng
include("./Coaraci/Coaraci.jl")

# Structure for plotting styles
struct MolSimStyle end
export MolSimStyle

end
