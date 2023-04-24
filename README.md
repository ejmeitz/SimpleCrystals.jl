# SimpleCrystals.jl

[![Build Status](https://ci.appveyor.com/api/projects/status/kd016pcm9epk1xk9?svg=true)](https://ci.appveyor.com/project/ejmeitz/simplecrystals-jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
<!-- [![Latest release](https://img.shields.io/github/release/ejmeitz/SimpleCrystals.jl.svg)](https://github.com/ejmeitz/SimpleCrystals.jl/releases/latest)
[![Documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMolSim.github.io/Molly.jl/stable)-->

 An interface for creating crystal geometries for molecular simulation. This package implements all monoatomic 3D and 2D Bravais lattices and enables to user to define a custom basis to create polyatomic Bravais lattice's or to create non-Bravais lattice structures.

 ### Examples

The functions to generate basic monoatomic crystal structures (i.e. fcc, bcc) are already implemented and require no customization.

##### Monoatomic Bravais Lattices:
Whenever possible crystals are implemented using a conventional unit cell so that patterning a simulation cell is simple. For FCC only the lattice parameter and element type are needed. SimpleCrystals.jl re-exports [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) and can handle any atomic symbols from [PeriodicTable.jl](https://github.com/JuliaPhysics/PeriodicTable.jl). The code below creates an FCC crystal of carbon with a conventional cell that is 5.4 Angstroms. The cell is patterened 4 times in the x, y, and z directions.

```julia
a = 0.54u"nm"
element = :C
fcc_crystal = FCC(a, element)
atoms = replicate_unit_cell(fcc_crystal, SVector(4,4,4))
```
![Monoatomic FCC](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png)

#### User Defined Crystal Structures
The SimpleCrystals API is not exhaustive, but provides an interface to create more complex, polyatomic crystals. For example NaCl




#### 3D Bravais Lattices
<table>
    <tr>
        <td>Crystal Family</td>
        <td>Primitive</td>
        <td>Base Centered</td>
        <td>Body Centered</td>
        <td>Face Centered</td>
    </tr>
    <tr>
        <td>Cubic</td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="1" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
    </tr>
    <tr>
        <td>Monoclinic</td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="1" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
    </tr>
    <tr>
        <td>Orthorhombic</td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="1" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
    </tr>
    <tr>
        <td>Tetragonal</td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="1" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
    </tr>
    <tr>
        <td>Hexagonal</td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="1" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
        <td> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 50px height = 32px> </td>
    </tr>

</table>

#### Other 3D Structrues
Diamond and HCP are also implemented as part of the API, but other non-bravais crystals can be created. 

BCC:
![BCC]((https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_bcc.png)
Diamond:
![Monoatomic Diamond](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_diamond.png)

Some crystals are not cubic and require a triclinic domain to properly replicate a unit cell, for example rhombohedral and HCP.

HCP:
![Monoatomic HCP](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_hcp.png)
Rhombohedral:
![Monoatomic Rhombohedral](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_rhomb.png)

#### 2D Bravais Lattices

#### Other 2D Structures
The honeycomb lattice is the only 2D non-bravais lattice implemented as part of the SimpleCrystals API.