module Reweighting

using CellListMap: ParticleSystem, map_pairwise, map_pairwise!
using ..MolSimToolkit: Simulation, positions, unitcell
using ..MolSimToolkit.MolecularMinimumDistances

export reweight, lennard_jones_perturbation, poly_decay_perturbation, gaussian_decay_perturbation

const testdir = "$(@__DIR__)/test"

include("./Reweight.jl")
include("./perturbation_examples.jl")

end #Module Reweighting

@testitem "Reweighting one frame" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_one_frame.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

    i2 = PDBTools.selindex(atoms(simulation), "residue 11")

    sum_of_dist = reweight(simulation, r -> r/10, [i1[239]], i2; cutoff = 25.0)
    @test sum_of_dist.energy ≈ [7.4295543149]
end

@testitem "Reweighting small trajectory" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "index 97 or index 106")

    i2 = PDBTools.selindex(atoms(simulation), "residue 15 and name HB3")

    sum_of_dist = reweight(simulation, r -> r/10, i1, i2, cutoff = 25.0)
    @test sum_of_dist.energy ≈ [
        1.773896547670759, 
        1.5923698293115915, 
        1.716614676290554, 
        1.933003841107648, 
        1.602329229247863, 
        1.9639005665480983, 
        3.573986006775934, 
        2.188798265022823, 
        2.066180657974777, 
        1.6845109623700647
    ]
end

@testitem "Reweighting small trajectory probs" begin
    import PDBTools
    using MolSimToolkit.Reweighting
    using MolSimToolkit.Reweighting: testdir

    simulation = Simulation("$testdir/Testing_reweighting.pdb", "$testdir/Testing_reweighting_10_frames_trajectory.xtc")

    i1 = PDBTools.selindex(atoms(simulation), "resname TFE and name O")

    i2 = PDBTools.selindex(atoms(simulation), "protein and name O")

    α = 5.e-3

    β = 5.e-3

    probs_test = reweight(simulation, r -> gaussian_decay_perturbation(r/10, α, β), i1, i2; cutoff = 12.0)
    @test probs_test.probability ≈ [
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
    ]
end