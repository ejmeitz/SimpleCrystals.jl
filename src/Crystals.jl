
export
    SC,
    FCC,
    BCC,
    replicate_unit_cell


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

#Convert to iterator

# Generates coordinates for 'crystal' from R = n1*a1 + n2*a2 + n3*a3 + basis
function replicate_unit_cell(crystal::Crystal{D}, N::SVector{D}) where D
    @assert all(N .> 0) "Number of unit cells should be positive"

    #Probably a way to get LP an not allocate all this memory
    # lattice_points = get_lattice_points(crystal.lattice, N)
    N_atoms = prod(N) * length(crystal.basis)

    #Create flat arrays for atoms & coords
    atoms = Vector{typeof(crystal.basis[1])}(undef,N_atoms)

    #Superimpose basis onto lattice points
    @inbounds for i in range(0,N_atoms-1)
        n = convert_1d_index(i, N)
        lattice_pt = crystal[n...]
        for basis_atom in crystal.basis
            atoms[i+1] = Atom(basis_atom.sym, lattice_pt .+ basis_atom.position, basis_atom.charge, basis_atom.mass)
        end
    end

    return atoms
end

#######################################################
### Monoatomic Examples with Conventional Unit Cell ###
#######################################################


# Monoatomic SC
function SC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a), zero(a), zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(a, zero(a), zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(zero(a), a, zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(zero(a), zero(a), a), charge = charge)]
    return Crystal(lattice,basis)
end

# Monoatomic FCC

function FCC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), FaceCentered())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(zero(a), 0.5*a, 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, zero(a), 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*a, zero(a)), charge = charge)]
    return Crystal(lattice, basis)
end

# Monoatomic BCC
function BCC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), BodyCentered())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis)
end


function Diamond(a, basis::Atom)

end