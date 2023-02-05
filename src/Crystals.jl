
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

#1D
function get_lattice_points(lattice::BravaisLattice{1}, N::SVector{1})

    Nx = N[1]
    lattice_points = SVector{Nx,SVector{D}}(undef)

    for i in range(0,Nx)
        lattice_points[i] = i.*lattice.primitive_vectors[1,:]
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
    coords = SVector{N_atoms,SVector{D}}

    #Superimpose basis onto lattice points
    i = 1
    for lp in lattice_points
        for basis_atom in crystal.basis
            coords[i] = lp .+ basis_atom.position
            atoms[i] = basis_atom.atom
            i += 1
        end
    end

    #Maximum extent of crystal in each direction
    dimensions = N .* crystal.lattice.lattice_constants




    return atoms, coords, dimensions
end


#Implement common crystal structures
struct FCC <: Crystal
    lattice::Lattice
    basis::SVector{BasisAtom}
end

#Monoatomic FCC
function FCC(a, atom::Atom)
    lattice = BravaisLattice(Cubic(a), FaceCentered())
    basis = SVector(BasisAtom(SVector(0.0,0.0,0.0), atom))
    return Crystal{3}(lattice,basis)
end

struct BCC <: Crystal
    
end

function BCC(a, atom::Atom)

end

struct Diamond <: Crystal

end

function Diamond(a, basis::Atom)

end