# SimpleCrystals.jl

[![Build Status](https://ci.appveyor.com/api/projects/status/kd016pcm9epk1xk9?svg=true)](https://ci.appveyor.com/project/ejmeitz/simplecrystals-jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
<!-- [![Latest release](https://img.shields.io/github/release/ejmeitz/SimpleCrystals.jl.svg)](https://github.com/ejmeitz/SimpleCrystals.jl/releases/latest)
[![Documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaMolSim.github.io/Molly.jl/stable)-->

 An interface for creating crystal geometries for molecular simulation. This package implements all monoatomic 3D and 2D Bravais lattices and enables to user to define a custom basis to create polyatomic Bravais lattices or to create non-Bravais crystal structures.

 ### Examples

The functions to generate basic monoatomic crystal structures (i.e. fcc, bcc) are already implemented and require no customization.

#### Monoatomic Bravais Lattices:
Whenever possible crystals are implemented using a conventional unit cell so that patterning a simulation cell is simple. A trinclinic boundary will work for the remaining lattices.

 For FCC only the lattice parameter and element type are needed. SimpleCrystals.jl re-exports [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) and can handle any atomic symbols from [PeriodicTable.jl](https://github.com/JuliaPhysics/PeriodicTable.jl). The code below creates an FCC crystal of carbon with a conventional cell that is 5.4 Angstroms. The cell is patterened 4 times in the x, y, and z directions.

```julia
a = 0.54u"nm"
element = :C
fcc_crystal = FCC(a, element)
atoms = get_coordinates(fcc_crystal, SVector(4,4,4))
```

#### User Defined Crystal Structures
The SimpleCrystals API is not exhaustive, but provides an interface to create non-bravais crystals and polyatomic crystals. For example, the Diamond crystal structure (which is a part of the API) is defined as simple cubic Bravais lattice with an 8 atom basis. Diamond is more naturally thought of as an FCC lattice with a 2 atom basis, but that would require a triclinic boundary.

The code below defines the BravaisLattice() object as a primitive, cubic lattice (simple cubic) with lattice parameter `a`. Then the basis is constructed as a list of Atom() objects. In this example, each basis atom is the same element but that could easily be changed. Finally, the Crystal() object is constructed from the BravaisLattice object and the list of basis atoms.

https://github.com/ejmeitz/SimpleCrystals.jl/blob/0ccc3f28e81d2c0aa5087039a52e94038520bad4/src/Crystals.jl#L99-L110

Similarly, we can create NaCl (not in the API) which can be thought of as two intertwined FCC lattices or a simple cubic lattice with a two atom basis.

Intertwined FCC:
```julia
a = 1.126u"nm"
```
2-Atom Basis SC:
```julia
a = 0.563u"nm" #half the lattice parameter of FCC example
```
Both methods yield the same list of coordinates
![NaCl Crystal]()


#### 3D Bravais Lattices
All 3D Bravais lattices created from the SimpleCrystal's API and visualized in [OVITO](https://ovito.org/).
<table>
    <tr>
        <th>Crystal Family</th>
        <th align="center">Primitive</th>
        <th align="center">Base Centered</th>
        <th align="center">Body Centered</th>
        <th align="center">Face Centered</th>
    </tr>
    <tr>
        <td align="center"><strong>Cubic</strong><br>a = b = c<br>&alpha; = &beta; = &gamma; = 90&#176</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/SC.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/BCC.png" alt="2" width = 160px height = 120px> </td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/FCC.png" alt="2" width = 160px height = 120px> </td>
    </tr>
    <tr>
        <td align="center"><strong>Triclinic</strong><br>a &#8800; b &#8800; c<br>&alpha; &#8800; &beta; &#8800; &gamma;</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/triclinic.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
        <td align="center"></td>
        <td align="center"></td>
    </tr>
    <tr>
        <td align="center"><strong>Monoclinic</strong><br>a &#8800; b &#8800; c<br>&alpha; = &gamma; = 90&#176, &beta;</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/monoclinic.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_base_centered.png" alt="2" width = 160px height = 120px> </td>
        <td align="center"></td>
        <td align="center"></td>
    </tr>
    <tr>
        <td align="center"><strong>Orthorhombic</strong><br>a &#8800; b &#8800; c<br>&alpha; = &beta; = &gamma; = 90&#176</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/ortho.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/ortho_base.png" alt="2" width = 160px height = 120px> </td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/ortho_body.png" alt="2" width = 160px height = 120px> </td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/ortho_fcc.png" alt="2" width = 160px height = 120px> </td>
    </tr>
    <tr>
        <td align="center"><strong>Tetragonal</strong><br>a = b &#8800; c<br> &alpha; = &beta; = &gamma; = 90&#176</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/tetragonal.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/tetragonal_body.png" alt="2" width = 160px height = 120px> </td>
        <td align="center"></td>
    </tr>
    <tr>
        <td align="center"><strong>Hexagonal (Rhombohedral)</strong><br>a = b = c<br>&alpha; = &beta; = &gamma;</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rhomb.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
        <td align="center"></td>
        <td align="center"></td>
    </tr>
      <tr>
        <td align="center"><strong>Hexagonal</strong><br>a = b &#8800; c<br>&alpha; = &beta; = 90&#176, &gamma; = 120&#176</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/hex_3d.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
        <td align="center"></td>
        <td align="center"></td>
    </tr>

</table>

#### Other 3D Structrues
Diamond and HCP are also implemented as part of the API: 
<table>
    <tr>
        <td align="center">Diamond</td>
        <td align="center">HCP</td>
    </tr>
    <tr>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_diamond.png" alt="2" width = 160px height = 120px> </td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_hcp.png" alt="2" width = 160px height = 120px> </td>
    </tr>
</table>


#### 2D Bravais Lattices
All 2D Bravais lattices created from the SimpleCrystal's API and visualized in [OVITO](https://ovito.org/).
<table>
    <tr>
        <th>Crystal Family</th>
        <th align="center">Primitive</th>
        <th align="center">Centered</th>
    </tr>
    <tr>
        <td >Monoclinic</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/oblique.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
    </tr>
    <tr>
        <td>Orthorhombic</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rect.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="2" width = 160px height = 120px> </td>
    </tr>
    <tr>
        <td>Tetragonal</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/square.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
    </tr>
    <tr>
        <td>Hexagonal</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/hex_2d.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
    </tr>

</table>


#### Other 2D Structures
The honeycomb lattice is the only 2D non-bravais lattice implemented as part of the SimpleCrystals API.
<table>
    <th align="center">Honeycomb</th>
    <tr>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_fcc.png" alt="1" width = 160px height = 120px> </td>
    </tr>
</table>

#### File I/O