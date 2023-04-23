
export
    SC,
    FCC,
    BCC,
    Diamond,
    HCP,
    Honeycomb,
    Rhombohedral,
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

# Generates coordinates for 'crystal' from R = n1*a1 + n2*a2 + n3*a3 + basis
function replicate_unit_cell(crystal::Crystal{D}, N::SVector{D}) where D
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

#######################################################
### Monoatomic Examples with Conventional Unit Cell ###
#######################################################
# All the strucrtures below use the conventional cell as this is the only way
# to pattern a cubic domain completely. This is consistent with the LAMMPS implementation.
# See https://github.com/lammps/lammps/blob/develop/src/lattice.cpp (85-116)

# Monoatomic SC -- 1 Atom Basis with SC Lattice Points
function SC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a), zero(a), zero(a)), charge = charge)]
    return Crystal(lattice,basis)
end

# Monoatomic FCC -- 4 Atom Basis with SC Lattice Points
function FCC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(0.5*a, 0.5*a, zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(zero(a), 0.5*a, 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, zero(a), 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice, basis)
end

# Monoatomic BCC -- 2 Atom Basis with SC Lattice Points
function BCC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*a, 0.5*a), charge = charge)]
    return Crystal(lattice,basis)
end

# Monoatomic Diamond -- 8 Atom Basis with SC Lattice Points
function Diamond(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(zero(a), 0.5*a, 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, zero(a), 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*a, zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.25*a, 0.25*a, 0.25*a), charge = charge),
             Atom(atomic_symbol, SVector(0.25*a, 0.75*a, 0.75*a), charge = charge),
             Atom(atomic_symbol, SVector(0.75*a, 0.25*a, 0.75*a), charge = charge),
             Atom(atomic_symbol, SVector(0.75*a, 0.75*a, 0.25*a), charge = charge)]
    return Crystal(lattice,basis)
end

###########################################################
# These crystals requrie a triclinc domain when building  #
# a simulation box for molecular dynamics                 #
###########################################################

# Hexagonal 2D Lattice -- 2 atom basis
function Honeycomb(a, atomic_symbol::Symbol; charge = 0.0u"C")
    d = a*sqrt(3)
    lattice = BravaisLattice(Hexagonal2D(a*sqrt(3)), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*d, 0.5*a), charge = charge),]
    return Crystal(lattice,basis)
end

# Hexagonal Lattice with a = c -- 2 atom basis
function HCP(a, atomic_symbol::Symbol; charge = 0.0u"C")
    c = a*sqrt(8/3)
    lattice = BravaisLattice(Hexagonal(a,c), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, a/(2*sqrt(3)), 0.5*c), charge = charge)]
    return Crystal(lattice,basis)
end

# Rhombohedral lattice -- 1 atom basis
# Simple cubic with α not equal to 90
function Rhombohedral(a, α, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Rhombohedral(a, α), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis)
end