var documenterSearchIndex = {"docs":
[{"location":"#SimpleCrystals.jl","page":"Home","title":"SimpleCrystals.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"(Image: Build Status) (Image: License: MIT) (Image: Latest release) (Image: Documentation stable)","category":"page"},{"location":"","page":"Home","title":"Home","text":"SimpleCrystals.jl is an interface for creating crystal geometries for molecular simulation within Julia. SimpleCrystals implements all 3D and 2D Bravais lattices (e.g. FCC & BCC) and allows users to define a custom basis to create polyatomic Bravais lattices or to create non-Bravais crystal structures. The AtomsBase interface is implemented to make use with other Julian molecular simulation software simple. There is no support for reading in crystal structures from other software (e.g. CIF files). SimpleCrystals is intended to provide a quick, lightweight method for generating atomic coordinates without leaving Julia.","category":"page"},{"location":"#Monoatomic-Bravais-Lattices:","page":"Home","title":"Monoatomic Bravais Lattices:","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Whenever possible crystals are implemented using a conventional unit cell so that patterning a simulation cell is simple. A trinclinic boundary will work for the remaining lattices.","category":"page"},{"location":"","page":"Home","title":"Home","text":"SimpleCrystals.jl re-exports Unitful.jl and StatcArrays.jl and can handle any atomic symbols from PeriodicTable.jl. For example, The code below creates an FCC crystal of carbon with a conventional cell that is 5.4 Angstroms. The cell is patterened 4 times in the x, y, and z directions. The coordinates can be obtained in code with the atoms member variable or saved to a XYZ file with to_xyz().","category":"page"},{"location":"","page":"Home","title":"Home","text":"a = 5.4u\"Å\"\nfcc_crystal = FCC(a, :C, SVector(4,4,4))\natoms = fcc_crystal.atoms\nto_xyz(fcc_crystal, raw\"./positions_fcc.xyz\")","category":"page"},{"location":"#D-Bravais-Lattices","page":"Home","title":"3D Bravais Lattices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"All 3D Bravais lattices created from the SimpleCrystal's API and visualized in OVITO. The radius of the atoms is chosen arbitrarily in OVITO.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<table>     <tr>         <th>Crystal Family</th>         <th align=\"center\">Primitive</th>         <th align=\"center\">Base Centered</th>         <th align=\"center\">Body Centered</th>         <th align=\"center\">Face Centered</th>     </tr>     <tr>         <td align=\"center\"><strong>Cubic</strong><br>a = b = c<br>&alpha; = &beta; = &gamma; = 90&#176</td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/SC.png\" alt=\"1\" width = 160px height = 120px> </td>         <td align=\"center\"></td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/BCC.png\" alt=\"2\" width = 160px height = 120px> </td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/FCC.png\" alt=\"2\" width = 160px height = 120px> </td>     </tr>     <tr>         <td align=\"center\"><strong>Triclinic</strong><br>a &#8800; b &#8800; c<br>&alpha; &#8800; &beta; &#8800; &gamma;</td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/triclinic.png\" alt=\"1\" width = 160px height = 120px> </td>         <td align=\"center\"></td>         <td align=\"center\"></td>         <td align=\"center\"></td>     </tr>     <tr>         <td align=\"center\"><strong>Monoclinic</strong><br>a &#8800; b &#8800; c<br>&alpha; = &gamma; = 90&#176, &beta;</td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/monoclinic.png\" alt=\"1\" width = 160px height = 120px> </td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/monobasecentered.png\" alt=\"2\" width = 160px height = 120px> </td>         <td align=\"center\"></td>         <td align=\"center\"></td>     </tr>     <tr>         <td align=\"center\"><strong>Orthorhombic</strong><br>a &#8800; b &#8800; c<br>&alpha; = &beta; = &gamma; = 90&#176</td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/ortho.png\" alt=\"1\" width = 160px height = 120px> </td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/orthobase.png\" alt=\"2\" width = 160px height = 120px> </td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/orthobody.png\" alt=\"2\" width = 160px height = 120px> </td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/orthofcc.png\" alt=\"2\" width = 160px height = 120px> </td>     </tr>     <tr>         <td align=\"center\"><strong>Tetragonal</strong><br>a = b &#8800; c<br> &alpha; = &beta; = &gamma; = 90&#176</td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/tetragonal.png\" alt=\"1\" width = 160px height = 120px> </td>         <td align=\"center\"></td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/tetragonalbody.png\" alt=\"2\" width = 160px height = 120px> </td>         <td align=\"center\"></td>     </tr>     <tr>         <td align=\"center\"><strong>Hexagonal (Rhombohedral)</strong><br>a = b = c<br>&alpha; = &beta; = &gamma;</td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rhomb.png\" alt=\"1\" width = 160px height = 120px> </td>         <td align=\"center\"></td>         <td align=\"center\"></td>         <td align=\"center\"></td>     </tr>       <tr>         <td align=\"center\"><strong>Hexagonal</strong><br>a, c<br>&alpha; = &beta; = 90&#176, &gamma; = 120&#176</td>         <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/hex_3d.png\" alt=\"1\" width = 160px height = 120px> </td>         <td align=\"center\"></td>         <td align=\"center\"></td>         <td align=\"center\"></td>     </tr>","category":"page"},{"location":"","page":"Home","title":"Home","text":"</table>","category":"page"},{"location":"#Other-3D-Structrues","page":"Home","title":"Other 3D Structrues","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Diamond and HCP are also implemented as part of the API: ","category":"page"},{"location":"","page":"Home","title":"Home","text":"<table>\n    <tr>\n        <td align=\"center\">Diamond</td>\n        <td align=\"center\">HCP</td>\n    </tr>\n    <tr>\n        <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/diamond.png\" alt=\"2\" width = 160px height = 120px> </td>\n        <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/HCP.png\" alt=\"2\" width = 160px height = 120px> </td>\n    </tr>\n</table>","category":"page"},{"location":"#D-Bravais-Lattices-2","page":"Home","title":"2D Bravais Lattices","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"All 2D Bravais lattices created from the SimpleCrystal's API and visualized in OVITO.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<table>\n    <tr>\n        <th>Crystal Family</th>\n        <th align=\"center\">Primitive</th>\n        <th align=\"center\">Centered</th>\n    </tr>\n    <tr>\n        <td><strong>Monoclinic</strong><br>a &#8800; b<br>&theta; &#8800; 90&#176</td>\n        <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/oblique.png\" alt=\"1\" width = 160px height = 120px> </td>\n        <td align=\"center\"></td>\n    </tr>\n    <tr>\n        <td><strong>Orthorhombic</strong><br>a &#8800; b<br>&theta; = 90&#176</td>\n        <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rect.png\" alt=\"1\" width = 160px height = 120px> </td>\n        <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/rect_centered.png\" alt=\"2\" width = 160px height = 120px> </td>\n    </tr>\n    <tr>\n        <td><strong>Tetragonal</strong><br>a = b<br>&theta; = 90&#176</td>\n        <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/square.png\" alt=\"1\" width = 160px height = 120px> </td>\n        <td align=\"center\"></td>\n    </tr>\n    <tr>\n        <td><strong>Hexagonal</strong><br>&theta; = 120&#176</td>\n        <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/hex_2d.png\" alt=\"1\" width = 160px height = 120px> </td>\n        <td align=\"center\"></td>\n    </tr>\n\n</table>","category":"page"},{"location":"#Other-2D-Structures","page":"Home","title":"Other 2D Structures","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The honeycomb lattice is the only 2D non-bravais lattice implemented as part of the SimpleCrystals API.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<table>\n    <th align=\"center\">Honeycomb</th>\n    <tr>\n        <td align=\"center\"> <img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/honeycomb.png\" alt=\"1\" width = 160px height = 120px> </td>\n    </tr>\n</table>","category":"page"},{"location":"#User-Defined-Crystal-Structures","page":"Home","title":"User Defined Crystal Structures","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The SimpleCrystals API is not exhaustive, but provides an interface to create non-bravais crystals and polyatomic crystals. For example, the Diamond crystal structure (which is a part of the API) is defined as simple cubic Bravais lattice with an 8 atom basis. Diamond is more naturally thought of as an FCC lattice with a 2 atom basis, but that would require a triclinic boundary.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The code below defines the BravaisLattice() object as a primitive, cubic lattice (simple cubic) with lattice parameter a. Then the basis is constructed as a list of Atom() objects. In this example, each basis atom is the same element but that could easily be changed. Finally, the Crystal() object is constructed from the BravaisLattice object and the list of basis atoms.","category":"page"},{"location":"","page":"Home","title":"Home","text":"https://github.com/ejmeitz/SimpleCrystals.jl/blob/0ccc3f28e81d2c0aa5087039a52e94038520bad4/src/Crystals.jl#L99-L110","category":"page"},{"location":"","page":"Home","title":"Home","text":"Similarly, we can create NaCl (not in the API) which can be thought of as a simple cubic lattice with an 8 atom basis or two intertwined FCC lattices.","category":"page"},{"location":"","page":"Home","title":"Home","text":"2-Atom Basis SC:","category":"page"},{"location":"","page":"Home","title":"Home","text":"function NaCl1(a, N::SVector{3})\n    lattice = BravaisLattice(Cubic(a), Primitive())\n    basis = [Atom(:Na, SVector(zero(a), zero(a), zero(a)), charge = 1.0u\"q\"),\n             Atom(:Na, SVector(0.5*a,zero(a),0.5*a), charge = 1.0u\"q\"),\n             Atom(:Na, SVector(zero(a), 0.5*a, 0.5*a), charge = 1.0u\"q\"),\n             Atom(:Na, SVector(0.5*a,0.5*a,zero(a)), charge = 1.0u\"q\"),\n             Atom(:Cl, SVector(0.5*a, zero(a), zero(a)), charge = -1.0u\"q\"),\n             Atom(:Cl, SVector(zero(a), 0.5*a, zero(a)), charge = -1.0u\"q\"),\n             Atom(:Cl, SVector(zero(a),zero(a),0.5*a), charge = -1.0u\"q\"),\n             Atom(:Cl, SVector(0.5*a, 0.5*a, 0.5*a), charge = -1.0u\"q\")]\n    return Crystal(lattice,basis,N)\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"Intertwined FCC:","category":"page"},{"location":"","page":"Home","title":"Home","text":"function NaCl2(a, N::SVector{3})\n    lattice = BravaisLattice(Cubic(a), FaceCentered())\n    basis = [Atom(:Na, SVector(zero(a), zero(a), zero(a)), charge = 1.0u\"q\"),\n             Atom(:Cl, SVector(0.5*a, zero(a), zero(a)), charge = -1.0u\"q\")]\n    return Crystal(lattice,basis,N)\nend","category":"page"},{"location":"","page":"Home","title":"Home","text":"Both methods yield the same structure with periodic boundary conditions, but the first function uses a conventional cell so the result is much easier to see and create a simulation box for. Whenever possible use a conventional cell (simple cubic lattice). Note that to use both of these functions the lattice parameter a is the distance between Na atoms (or Cl atoms) not the Na-Cl distance as the basis places the atoms at the proper 0.5*a spacing.","category":"page"},{"location":"","page":"Home","title":"Home","text":"<table>\n<tr>\n    <th align=\"center\">Conventional Cell</th>\n    <th align=\"center\">FCC Unit Cell</th>\n</tr>\n<tr>\n    <td align=\"center\"><img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/NaCl_8atom_basis.png\" alt=\"1\" width = 320px height = 240px></td>\n    <td align=\"center\"><img src=\"https://github.com/ejmeitz/SimpleCrystals.jl/raw/main/assets/nacl_fcc_basis.png\" alt=\"1\" width = 320px height = 240px></td>\n</tr>\n</table>","category":"page"},{"location":"#File-I/O","page":"Home","title":"File I/O","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Using one of the built-in crystal objects (e.g. FCC) or a user-defined crystal you can call the to_xyz() function to print out a .xyz file with a chosen number of unit cells. For example, to get the coordinates for 4 unit-cells in the x, y and z directions for FCC you could use the following code:","category":"page"},{"location":"","page":"Home","title":"Home","text":"a = 5.4u\"Å\"\nfcc_crystal = FCC(a, :C, SVector(4,4,4))\nto_xyz(fcc_crystal, raw\"C:\\Users\\ejmei\\Desktop\\positions_fcc.xyz\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"To get the list of atoms in code you can use the atoms member variable. The code below will return the array of atoms that the to_xyz() function builds internally.","category":"page"},{"location":"","page":"Home","title":"Home","text":"a = 5.4u\"Å\"\nfcc_crystal = FCC(a, :C, SVector(4,4,4))\natoms = fcc_crystal.atoms","category":"page"}]
}
