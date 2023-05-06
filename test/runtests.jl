using SimpleCrystals
using Test
using LinearAlgebra

# These tests are only designed to ensure that the constructors for
# the crystals work. It does not check that all the coordinates
# generated are correct.

@testset "3D-Bravais" begin
    a = 0.54u"nm"
    sc_crystal = SC(0.54u"nm", :C, SVector(4,4,4))
    fcc_crystal = FCC(0.54u"nm", :C, SVector(4,4,4))
    diamond_crystal = Diamond(a, :C, SVector(4,4,4))
    hex_crystal = Honeycomb(0.54u"nm", :C, SVector(4,4))


    @test length(atoms_sc) == 64
    @test norm(ustrip.(sc_crystal[1].position) .- ustrip.(sc_crystal[2].position)) == ustrip(a)

    @test length(atoms_fcc) == 256
    @test norm(ustrip.(fcc_crystal[1].position) .- ustrip.(fcc_crystal[2].position)) ≈ ustrip(a)/sqrt(2)

    @test length(diamond_crystal) == 512

    @test length(hex_crystal) == 32

    # Add tests to check that primitive vectors are right handed
end

@testset "AtomsBase" begin
    
    fcc_crystal = FCC(0.54u"nm", :C, SVector(4,4,4))

    @test n_dimensions(fcc_crystal) == 3
    @test atomic_mass(fcc_crystal, 1) == atomic_mass(fcc_crystal[1])
    @test atomic_number(fcc_crystal,1) == :C
end


