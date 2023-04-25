using SimpleCrystals
using Test
using LinearAlgebra

# These tests are only designed to ensure that the constructors for
# the crystals work. It does not check that all the coordinates
# generated are correct.

@testset "3D-Bravais" begin
    a = 0.54u"nm"
    sc_crystal = SC(0.54u"nm", :C)
    atoms_sc = get_coordinates(sc_crystal, SVector(4,4,4))

    fcc_crystal = FCC(0.54u"nm", :C)
    atoms_fcc = get_coordinates(fcc_crystal, SVector(4,4,4))

    @test length(atoms_sc) == 64
    @test norm(ustrip.(atoms_sc[1].position) .- ustrip.(atoms_sc[2].position)) == ustrip(a)

    @test length(atoms_fcc) == 256
    @test norm(ustrip.(atoms_fcc[1].position) .- ustrip.(atoms_fcc[2].position)) ≈ ustrip(a)/sqrt(2)
end

@testset "3D-Other" begin
    

end

@testset "2D-Bravais" begin
    

end

@testset "2D-Other" begin
    
    hex_crystal = Honeycomb(0.54u"nm",:C)
    atoms = get_coordinates(hex_crystal, SVector(4,4))

    @test length(atoms) == 32

end


