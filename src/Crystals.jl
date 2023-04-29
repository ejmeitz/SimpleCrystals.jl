
using AtomsBase

export Crystal

struct Crystal{D, A, B <: AbstractVector{<:Atom{D}}} #<: AbstractSystem{D}
    lattice::BravaisLattice{D}
    basis::B
    N_unit_cells::SVector{D}
    atoms::A
end

#The type inference for D is kind of jank
function Crystal(lattice, basis, N_unit_cells)
    atoms = get_coordinates(lattice, basis, N_unit_cells)
    return Crystal{size(lattice.primitive_vectors)[1], typeof(basis), typeof(atoms)}(lattice, basis, N_unit_cells, atoms)
end

function num_atoms_per_unit_cell(crystal::Crystal)
    return length(crystal.basis)
end


#Returns R = n1*a1 + n2*a2 + n3*a3, where a are the primitive lattice vectors
function Base.getindex(lattice::BravaisLattice{D}, indices::Vararg{Integer,D}) where D
    return SVector{D}(sum(indices .* lattice.primitive_vectors, dims = 1))
end

function Base.getindex(crystal::Crystal, i::Int)
    return crystal.atoms[i]
end

# Convert 1D index to equivalent 3D index in a matrix with dimensions N
# Note this is 0-indexed
function convert_1d_index(i, N::SVector{3})
    z = i ÷ (N[1] * N[2])
    i -= (z * N[1] * N[2])
    y = i ÷ N[1]
    x = i % N[1]
    return x, y, z
end

# Convert 1D index to equivalent 2D index in a matrix with dimensions N
# Note this is 0-indexed
function convert_1d_index(i, N::SVector{2})
    x = i % N[1]
    y = i ÷ N[1]
    return x, y
end

# Generates coordinates for 'crystal' from R = n1*a1 + n2*a2 + n3*a3 + basis
function get_coordinates(lattice::BravaisLattice{D}, basis, N::SVector{D}) where D
    @assert all(N .> 0) "Number of unit cells should be greater than 0"

    #Probably a way to get LP an not allocate all this memory
    N_lattice_pts = prod(N)
    N_basis_atoms = length(basis)
    N_atoms = N_lattice_pts * N_basis_atoms

    #Create flat arrays for atoms & coords
    atoms = Vector{typeof(basis[1])}(undef,N_atoms)

    #Superimpose basis onto lattice points
    @inbounds for i in range(0,N_lattice_pts-1)
        n = convert_1d_index(i, N)
        lattice_pt = lattice[n...]
        for (j,basis_atom) in enumerate(basis)
            atoms[N_basis_atoms*i + j] = Atom(basis_atom.sym, lattice_pt .+ basis_atom.position, basis_atom.charge, basis_atom.mass)
        end
    end

    return atoms
end


### AtomsBase Compliance ###

Base.length(sys::Crystal) = length(sys.atoms)
# AtomsBase.bounding_box(sys::Crystal) = sys.system_data.bounding_box
# AtomsBase.boundary_conditions(sys::Crystal) = sys.system_data.boundary_conditions

# AtomsBase.atomkeys(sys::Crystal) = keys(sys.atom_data)
# AtomsBase.hasatomkey(sys::Crystal, x::Symbol) = haskey(sys.atom_data, x)

# AtomsBase.position(sys::Crystal)                  = 
AtomsBase.position(sys::Crystal, i::Integer)        = sys.atoms[i].position
# AtomsBase.velocity(sys::Crystal)                  = Base.getindex(s, :, :velocity)
# AtomsBase.velocity(sys::Crystal, i::Integer)      = Base.getindex(s, i, :velocity)
# AtomsBase.atomic_mass(sys::Crystal)               = Base.getindex(s, :, :atomic_mass)
# AtomsBase.atomic_mass(sys::Crystal, i::Integer)   = Base.getindex(s, i, :atomic_mass)
# AtomsBase.atomic_symbol(sys::Crystal)             = Base.getindex(s, :, :atomic_symbol)
# AtomsBase.atomic_symbol(sys::Crystal, i::Integer) = Base.getindex(s, i, :atomic_symbol)
# AtomsBase.atomic_number(sys::Crystal)             = Base.getindex(s, :, :atomic_number)
# AtomsBase.atomic_number(sys::Crystal, i::Integer) = Base.getindex(s, i, :atomic_number)