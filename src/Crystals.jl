
using AtomsBase

export Crystal

struct Crystal{D, A, B <: AbstractVector{<:Atom{D}}} <: AbstractSystem{D}
    lattice::BravaisLattice{D}
    basis::B
    N_unit_cells::SVector{D}
    atoms::A
end

#The type inference for D is kind of jank
# Override of Base functions at bottom of file
function Crystal(lattice, basis, N_unit_cells)
    atoms = get_coordinates(lattice, basis, N_unit_cells)
    return Crystal{size(lattice.primitive_vectors)[1], typeof(basis), typeof(atoms)}(lattice, basis, N_unit_cells, atoms)
end


function atoms_per_unit_cell(crystal::Crystal)
    return length(crystal.basis)
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

# Generates coordinates for a crystal struct from R = n1*a1 + n2*a2 + n3*a3 + basis
function get_coordinates(lattice::BravaisLattice{D}, basis, N::SVector{D}) where D
    @assert all(N .> 0) "Number of unit cells should be greater than 0"

    #Probably a way to get LP an not allocate all this memory
    N_lattice_pts = prod(N)
    N_basis_atoms = length(basis)
    N_atoms = N_lattice_pts * N_basis_atoms

    #Create flat arrays for atoms & coords
    atoms = Vector{typeof(basis[1])}(undef,N_atoms)

    #Superimpose basis onto lattice points
    for i in range(0,N_lattice_pts-1)
        n = convert_1d_index(i, N)
        lattice_pt = lattice[n...]
        for (j,basis_atom) in enumerate(basis)
            atoms[N_basis_atoms*i + j] = Atom(basis_atom.sym, lattice_pt .+ basis_atom.position, basis_atom.charge, basis_atom.mass)
        end
    end

    return atoms
end


### AtomsBase Compliance ###

Base.getindex(sys::Crystal, i::Int) = sys.atoms[i]
Base.getindex(sys::Crystal, i::Int, x::Symbol) = Base.getindex(sys.atoms[i], x)
Base.getindex(sys::Crystal, ::Colon, x::Symbol) = Base.getindex.(sys.atoms, Ref(x))
AtomsBase.keys(sys::Crystal) = (:lattice, :basis, :N_unit_cells, :atoms)
Base.getindex(sys::Crystal, x::Symbol) = hasfield(Crystal, x) ? getfield(sys, x) : error("No field `$x` in Atom object. Allowed keys are $(keys(sys)).")
AtomsBase.haskey(sys::Crystal, x::Symbol) = hasfield(Crystal, x)
Base.pairs(sys::Crystal) = (k => sys[k] for k in keys(sys))


Base.length(sys::Crystal) = length(sys.atoms)
Base.iterate(sys::Crystal, state = 1) = state > length(sys) ? nothing : (sys.atoms[state], state + 1)
Base.eachindex(sys::Crystal) = Base.OneTo(length(sys))

AtomsBase.bounding_box(sys::Crystal) = sys.lattice.primitive_vectors .* sys.N_unit_cells
AtomsBase.boundary_conditions(sys::Crystal{3}) = SVector(Periodic(),Periodic(),Periodic())
AtomsBase.boundary_conditions(sys::Crystal{2}) = SVector(Periodic(),Periodic())
AtomsBase.n_dimensions(sys::Crystal) = length(sys.N_unit_cells)

AtomsBase.atomkeys(sys::Crystal) = keys(sys.atoms[1])
AtomsBase.hasatomkey(sys::Crystal, x::Symbol) = haskey(sys.atoms[1], x)

AtomsBase.species_type(sys::Crystal) = typeof(sys.atoms[1])

AtomsBase.position(sys::Crystal) = position.(sys.atoms)
AtomsBase.position(sys::Crystal, i::Integer) = sys.atoms[i].position
AtomsBase.position(sys::Crystal, ::Colon) = position.(sys.atoms)

AtomsBase.velocity(sys::Crystal) = missing #zeros(size(sys.atoms)) * u"m * s^-1"
AtomsBase.velocity(sys::Crystal, i::Integer) = missing #zeros(size(sys.atoms[i]))* u"m * s^-1"
AtomsBase.velocity(sys::Crystal, ::Colon) = missing

AtomsBase.atomic_mass(sys::Crystal) = atomic_mass.(sys.atoms)
AtomsBase.atomic_mass(sys::Crystal, i::Integer) = sys.atoms[i].mass
AtomsBase.atomic_mass(sys::Crystal, ::Colon) = atomic_mass.(sys.atoms)

AtomsBase.atomic_symbol(sys::Crystal) = atomic_symbol.(sys.atoms)
AtomsBase.atomic_symbol(sys::Crystal, i::Integer) = sys.atoms[i].sym

AtomsBase.atomic_number(sys::Crystal, ::Colon) = atomic_number.(sys.atoms)
AtomsBase.atomic_number(sys::Crystal) = atomic_number.(sys.atoms)
AtomsBase.atomic_number(sys::Crystal, i::Integer) = atomic_number(sys.atoms[i])

AtomsBase.visualize_ascii(sys::Crystal) = ""
Base.show(io::IO, sys::Crystal) = print(io, "Crystal with ", length(sys)," atoms.")