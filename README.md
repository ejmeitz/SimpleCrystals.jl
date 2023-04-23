# SimpleCrystals.jl

<!-- [![CI](https://ci.appveyor.com/project/ejmeitz/SimpleCrystals.jl/branch/master)]()
[![Latest release](https://img.shields.io/github/release/ejmeitz/SimpleCrystals.jl.svg)](https://github.com/ejmeitz/SimpleCrystals.jl/releases/latest) -->
[![License](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/ejmeitz/SimpleCrystals.jl/blob/master/LICENSE.md)

 An interface for creating crystal geometries for molecular simulation. This package is not intended to provide a comprehensive system for generating all types of crstyals. The main goal is to enable users to rapidly generate atomic coordinates for common crystal structures like FCC and BCC without much effort.

 ### Examples

The functions to generate basic monoatomic crystal structures (i.e. fcc, bcc) are already implemented and require no customization.

##### Monoatomic FCC
Whenever possible crystals are implemented using a conventional unit cell so that patterning a simulation cell is simple. For FCC only the lattice parameter and element type are needed. SimpleCrystals.jl re-exports [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) and can handle any atomic symbols from [PeriodicTable.jl](https://github.com/JuliaPhysics/PeriodicTable.jl). The code below creates an FCC crystal of Carbon with a conventional cell that is 5.4 Angstroms. The cell is patterened 4 times in the x, y, and z directions.

'''
a = 0.54u"nm"
element = :C
fcc_crystal = FCC(0.54u"nm", element)
atoms = replicate_unit_cell(fcc_crystal, SVector(4,4,4))
'''
