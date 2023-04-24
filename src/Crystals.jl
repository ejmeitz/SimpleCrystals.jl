
export
    get_coordinates


struct Crystal{D, B <: AbstractVector{<:Atom{D}}}
    lattice::BravaisLattice{D}
    basis::B
end

#Returns R = n1*a1 + n2*a2 + n3*a3, where a are the primitive lattice vectors
    #This will just be a lattice point, does not account for basis
function Base.getindex(crystal::Crystal{D}, indices::Vararg{Integer,D}) where D
    return SVector{D}(sum(indices .* crystal.lattice.primitive_vectors, dims = 1))
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
function get_coordinates(crystal::Crystal{D}, N::SVector{D}) where D
    @assert all(N .> 0) "Number of unit cells should be greater than 0"

    #Probably a way to get LP an not allocate all this memory
    N_lattice_pts = prod(N)
    N_basis_atoms = length(crystal.basis)
    N_atoms = N_lattice_pts * N_basis_atoms

    #Create flat arrays for atoms & coords
    atoms = Vector{typeof(crystal.basis[1])}(undef,N_atoms)

    #Superimpose basis onto lattice points
    @inbounds for i in range(0,N_lattice_pts-1)
        n = convert_1d_index(i, N)
        lattice_pt = crystal[n...]
        for (j,basis_atom) in enumerate(crystal.basis)
            atoms[N_basis_atoms*i + j] = Atom(basis_atom.sym, lattice_pt .+ basis_atom.position, basis_atom.charge, basis_atom.mass)
        end
    end

    return atoms
end
