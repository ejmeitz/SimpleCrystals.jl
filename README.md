# SimpleCrystals.jl

[![Build Status](https://ci.appveyor.com/api/projects/status/kd016pcm9epk1xk9?svg=true)](https://ci.appveyor.com/project/ejmeitz/simplecrystals-jl)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![Latest release](https://img.shields.io/github/release/ejmeitz/SimpleCrystals.jl.svg)](https://github.com/ejmeitz/SimpleCrystals.jl/releases/latest)
[![Documentation stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ejmeitz.github.io/SimpleCrystals.jl/stable/)

 SimpleCrystals.jl is an interface for creating crystal geometries for molecular simulation within Julia. SimpleCrystals implements all 3D and 2D Bravais lattices (e.g. FCC & BCC) and allows users to define a custom basis to create polyatomic Bravais lattices or to create non-Bravais crystal structures. The [AtomsBase](https://github.com/JuliaMolSim/AtomsBase.jl) interface is implemented to make use with other Julian molecular simulation software simple. There is no support for reading in crystal structures from other software (e.g. CIF files). Check out [Xtals.jl](https://github.com/SimonEnsemble/Xtals.jl) or [AtomsIO.jl](https://github.com/mfherbst/AtomsIO.jl) for that. SimpleCrystals is intended to provide a quick, lightweight method for generating atomic coordinates without leaving Julia.



#### Monoatomic Bravais Lattices:
Whenever possible crystals are implemented using a conventional unit cell so that patterning a simulation cell is simple. A trinclinic boundary will work for the remaining lattices.

 SimpleCrystals.jl re-exports [Unitful.jl](https://painterqubits.github.io/Unitful.jl/stable/) and [StatcArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) and can handle any atomic symbols from [PeriodicTable.jl](https://github.com/JuliaPhysics/PeriodicTable.jl). For example, the code below creates an FCC crystal of carbon with a conventional cell that is 5.4 Angstroms. The cell is patterened 4 times in the x, y, and z directions. The coordinates can be obtained in code with the `atoms` member variable or saved to a XYZ file with `to_xyz()`.

```julia
a = 5.4u"Å"
fcc_crystal = FCC(a, :C, SVector(4,4,4))
atoms = fcc_crystal.atoms
to_xyz(fcc_crystal, raw"./positions_fcc.xyz")

#Equivalently if you do not want to specify an atomic species
a = 5.4u"Å"
fcc_crystal = FCC(a, 12.01u"g/mol", SVector(4,4,4))
atoms = fcc_crystal.atoms
to_xyz(fcc_crystal, raw"./positions_fcc.xyz")
```

#### 3D Bravais Lattices
All 3D Bravais lattices created from the SimpleCrystal's API and visualized in [OVITO](https://ovito.org/). The radius of the atoms is chosen arbitrarily in OVITO.
The full list of implemented functions can be found [here](https://github.com/ejmeitz/SimpleCrystals.jl/blob/main/src/bravais/3D_bravais.jl). 

| Crystal Family | Primitive | Base Centered | Body Centered | Face Centered |
|     :---:      |   :---:   |     :---:     |     :---:     |     :---:     |
| **Cubic**<br>a = b = c<br>α = β = γ = 90° | ![1](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/SC.png) |  | ![2](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/BCC.png) | ![2](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/FCC.png) |
| **Triclinic**<br>a ≠ b ≠ c<br>α ≠ β ≠ γ | ![1](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/triclinic.png) |  |  |  |
| **Monoclinic**<br>a ≠ b ≠ c<br>α = γ = 90°, β | ![1](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/monoclinic.png) | ![2](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/mono_base_centered.png) |  |  |
| **Orthorhombic**<br>a ≠ b ≠ c<br>α = β = γ = 90° | ![1](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/ortho.png) | ![2](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/ortho_base.png) | ![2](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/ortho_body.png) | ![2](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/ortho_fcc.png) |
| **Tetragonal**<br>a = b ≠ c<br>α = β = γ = 90° | ![1](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/tetragonal.png) |  | ![2](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/tetragonal_body.png) |  |
| **Hexagonal (Rhombohedral)**<br>a = b = c<br>α = β = γ | ![1](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rhomb.png) |  |  |  |
| **Hexagonal**<br>a, c<br>α = β = 90°, γ = 120° | ![1](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/hex_3d.png) |  |  |  |






<!-- <table>
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
        <td align="center"><strong>Hexagonal</strong><br>a, c<br>&alpha; = &beta; = 90&#176, &gamma; = 120&#176</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/hex_3d.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
        <td align="center"></td>
        <td align="center"></td>
    </tr>

</table> -->

#### Other 3D Structrues
Diamond and HCP are also implemented as part of the API: 

| Diamond | HCP |
|:-------:|:---:|
| ![Diamond](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/diamond.png) | ![HCP](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/HCP.png) |


<!-- <table>
    <tr>
        <td align="center">Diamond</td>
        <td align="center">HCP</td>
    </tr>
    <tr>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/diamond.png" alt="2" width = 160px height = 120px> </td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/HCP.png" alt="2" width = 160px height = 120px> </td>
    </tr>
</table> -->


#### 2D Bravais Lattices
All 2D Bravais lattices created from the SimpleCrystal's API and visualized in [OVITO](https://ovito.org/).
The full list of implemented functions can be found [here](https://github.com/ejmeitz/SimpleCrystals.jl/blob/main/src/bravais/2D_bravais.jl). 

| Crystal Family | Primitive | Centered |
|:--------------:|:---------:|:--------:|
| **Monoclinic**<br>a ≠ b<br>θ ≠ 90° | ![Primitive](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/oblique.png) |  |
| **Orthorhombic**<br>a ≠ b<br>θ = 90° | ![Primitive](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rect.png) | ![Centered](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rect_centered.png) |
| **Tetragonal**<br>a = b<br>θ = 90° | ![Primitive](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/square.png) |  |
| **Hexagonal**<br>θ = 120° | ![Primitive](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/hex_2d.png) |  |

<!-- 
<table>
    <tr>
        <th>Crystal Family</th>
        <th align="center">Primitive</th>
        <th align="center">Centered</th>
    </tr>
    <tr>
        <td><strong>Monoclinic</strong><br>a &#8800; b<br>&theta; &#8800; 90&#176</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/oblique.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
    </tr>
    <tr>
        <td><strong>Orthorhombic</strong><br>a &#8800; b<br>&theta; = 90&#176</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rect.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rect_centered.png" alt="2" width = 160px height = 120px> </td>
    </tr>
    <tr>
        <td><strong>Tetragonal</strong><br>a = b<br>&theta; = 90&#176</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/square.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
    </tr>
    <tr>
        <td><strong>Hexagonal</strong><br>&theta; = 120&#176</td>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/hex_2d.png" alt="1" width = 160px height = 120px> </td>
        <td align="center"></td>
    </tr>

</table> -->


#### Other 2D Structures
The honeycomb lattice is the only 2D non-bravais lattice implemented as part of the SimpleCrystals API.

| Honeycomb |
|:---------:|
| ![Honeycomb](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/honeycomb.png) |


<!-- <table>
    <th align="center">Honeycomb</th>
    <tr>
        <td align="center"> <img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/honeycomb.png" alt="1" width = 160px height = 120px> </td>
    </tr>
</table> -->


#### User Defined Crystal Structures
The SimpleCrystals API is not exhaustive, but provides an interface to create non-bravais crystals and polyatomic crystals. For example, the Diamond crystal structure (which is a part of the API) is defined as simple cubic Bravais lattice with an 8 atom basis. Diamond is more naturally thought of as an FCC lattice with a 2 atom basis, but that would require a triclinic boundary.

The code below defines the BravaisLattice() object as a primitive, cubic lattice (simple cubic) with lattice parameter `a`. Then the basis is constructed as a list of Atom() objects. In this example, each basis atom is the same element but that could easily be changed. Finally, the Crystal() object is constructed from the BravaisLattice object and the list of basis atoms.

```julia
function Diamond(a, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)],atomic_mass, charge = charge),
             Atom([zero(a), 0.5*a, 0.5*a], atomic_mass, charge = charge),
             Atom([0.5*a, zero(a), 0.5*a], atomic_mass, charge = charge),
             Atom([0.5*a, 0.5*a, zero(a)], atomic_mass, charge = charge),
             Atom([0.25*a, 0.25*a, 0.25*a], atomic_mass, charge = charge),
             Atom([0.25*a, 0.75*a, 0.75*a], atomic_mass, charge = charge),
             Atom([0.75*a, 0.25*a, 0.75*a], atomic_mass, charge = charge),
             Atom([0.75*a, 0.75*a, 0.25*a], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end
```
Similarly, we can create NaCl (not in the API) which can be thought of as a simple cubic lattice with an 8 atom basis or two intertwined FCC lattices.

2-Atom Basis SC:
```julia
function NaCl1(a, N::SVector{3})
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom(:Na, [zero(a), zero(a), zero(a)], charge = 1.0u"q"),
             Atom(:Na, [0.5*a,zero(a),0.5*a], charge = 1.0u"q"),
             Atom(:Na, [zero(a), 0.5*a, 0.5*a], charge = 1.0u"q"),
             Atom(:Na, [0.5*a,0.5*a,zero(a)], charge = 1.0u"q"),
             Atom(:Cl, [0.5*a, zero(a), zero(a)], charge = -1.0u"q"),
             Atom(:Cl, [zero(a), 0.5*a, zero(a)], charge = -1.0u"q"),
             Atom(:Cl, [zero(a),zero(a),0.5*a], charge = -1.0u"q"),
             Atom(:Cl, [0.5*a, 0.5*a, 0.5*a], charge = -1.0u"q")]
    return Crystal(lattice,basis,N)
end
```
Intertwined FCC:
```julia
function NaCl2(a, N::SVector{3})
    lattice = BravaisLattice(CubicLattice(a), FaceCentered())
    basis = [Atom(:Na, [zero(a), zero(a), zero(a)], charge = 1.0u"q"),
             Atom(:Cl, [0.5*a, zero(a), zero(a)], charge = -1.0u"q")]
    return Crystal(lattice,basis,N)
end
```

Both methods yield the same structure with periodic boundary conditions, but the first function uses a conventional cell so the result is much easier to see and create a simulation box for. Whenever possible use a conventional cell (simple cubic lattice). Note that to use both of these functions the lattice parameter a is the distance between Na atoms (or Cl atoms) not the Na-Cl distance as the basis places the atoms at the proper 0.5*a spacing.

| Conventional Cell | FCC Unit Cell |
|:-----------------:|:--------------:|
| ![Conventional Cell](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/NaCl_8atom_basis.png) | ![FCC Unit Cell](https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/nacl_fcc_basis.png) |


<!-- <table>
<tr>
    <th align="center">Conventional Cell</th>
    <th align="center">FCC Unit Cell</th>
</tr>
<tr>
    <td align="center"><img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/NaCl_8atom_basis.png" alt="1" width = 320px height = 240px></td>
    <td align="center"><img src="https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/nacl_fcc_basis.png" alt="1" width = 320px height = 240px></td>
</tr>
</table> -->

#### File I/O

Using one of the built-in crystal objects (e.g. FCC) or a user-defined crystal you can call the `to_xyz()` function to print out a .xyz file with a chosen number of unit cells. For example, to get the coordinates for 4 unit-cells in the x, y and z directions for FCC you could use the following code:

```julia
a = 5.4u"Å"
fcc_crystal = FCC(a, :C, SVector(4,4,4))
to_xyz(fcc_crystal, raw"C:\Users\ejmei\Desktop\positions_fcc.xyz")
```

To get the list of atoms in code you can use the `atoms` member variable. The code below will return the array of atoms that the to_xyz() function builds internally. Alternatively, you can loop over the crystal object which will iterate through the `crystal.atoms` list.

```julia
a = 5.4u"Å"
fcc_crystal = FCC(a, :C, SVector(2,2,2))
atoms = fcc_crystal.atoms

for atom in fcc_crystal
    #do something
end
```
