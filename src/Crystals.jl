
export
    FCC,
    BCC


struct Crystal{D,L}
    lattice::BravaisLattice{D}
    basis::SVector{Atom}
    length_unit::L
end

function Crystal(lattice, basis)

end

# 3D version
function get_lattice_points(lattice::BravaisLattice{3}, N::SVector{3})

    lattice_points = SVector{prod(N),SVector{D}}(undef)
    Nx = N[1]
    Ny = N[2]
    Nz = N[3]

    for i in range(0, Nx), j in range(0, Ny), k in range(0, Nz)
        idx = (i + (j * Nx) + (k * Nx * Ny)) + 1 #plus 1 because of 1-indexed arrays
        lattice_points[idx] = i.*lattice.primitive_vectors[1,:] .+ j.*lattice.primitive_vectors[2,:] .+ k.lattice.primitive_vectors[3,:]
    end

    return lattice_points
end

#2D version
function get_lattice_points(lattice::BravaisLattice{2}, N::SVector{2})


    lattice_points = SVector{prod(N),SVector{D}}(undef)
    Nx = N[1]
    Ny = N[2]

    for i in range(0,Nx), j in range(0,Ny)
        idx = (i + (j * Nx)) + 1 #plus 1 because of 1-indexed arrays
        lattice_points[idx] = i.*lattice.primitive_vectors[1,:] .+ j.*lattice.primitive_vectors[2,:]
    end

    return lattice_points
end


function replicate_unit_cell(crystal::Crystal{D}, N::SVector{D}) where D
    @assert all(N .> 0) "Number of unit cells should be positive"

    #Probably a way to get LP an not allocate all this memory
    lattice_points = get_lattice_points(crystal.lattice, N)
    N_atoms = length(lattice_points) * length(crystal.basis)

    #Create flat arrays for atoms & coords
    atoms = SVector{N_atoms,Atom}

    #Superimpose basis onto lattice points
    i = 1
    for lp in lattice_points
        for basis_atom in crystal.basis         
            atoms[i] = Atom(basis_atom, lp .+ basis_atom.position)
            i += 1
        end
    end

    return atoms
end


# #Implement common crystal structures
# struct FCC <: Crystal
#     lattice::BravaisLattice{D}
#     basis::SVector{BasisAtom{D}}
# end

#Monoatomic FCC

function FCC(a, atom::Atom)
    lattice = BravaisLattice(Cubic(a), FaceCentered())
    basis = SVector(Atom(atom, SVector(zero(a),zero(a),zero(a))))
    return Crystal{3}(lattice, basis)
end

# struct BCC <: Crystal
    
# end

function BCC(a, atom::Atom)
    lattice = BravaisLattice(Cubic(a), BodyCentered())
    basis = SVector(atom)
    return Crystal{3}(lattice,basis)
end

# struct Diamond <: Crystal

# end

function Diamond(a, basis::Atom)

end