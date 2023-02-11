
export
    FCC,
    BCC


struct Crystal{D, B <: AbstractVector{<:Atom{D}}}
    lattice::BravaisLattice{D}
    basis::B
end

#Returns R = n1*a1 + n2*a2 + n3*a3, where a are the primitive lattice vectors
    #This will just be a lattice point, does not account for basis
function Base.getindex(crystal::Crystal{D}, indices::Vararg{Integer,D})
    return sum(indices .* crystal.lattice.primitive_vectors, dims = 1)
end



# # 3D version
# function get_lattice_points(lattice::BravaisLattice{3}, N::SVector{3})

#     lattice_points = SVector{prod(N),SVector{D}}(undef)
#     Nx = N[1]
#     Ny = N[2]
#     Nz = N[3]

#     for i in range(0, Nx), j in range(0, Ny), k in range(0, Nz)
#         idx = (i + (j * Nx) + (k * Nx * Ny)) + 1 #plus 1 because of 1-indexed arrays
#         lattice_points[idx] = i.*lattice.primitive_vectors[1,:] .+ j.*lattice.primitive_vectors[2,:] .+ k.lattice.primitive_vectors[3,:]
#     end

#     return lattice_points
# end

# #2D version
# function get_lattice_points(lattice::BravaisLattice{2}, N::SVector{2})


#     lattice_points = SVector{prod(N),SVector{D}}(undef)
#     Nx = N[1]
#     Ny = N[2]

#     for i in range(0,Nx), j in range(0,Ny)
#         idx = (i + (j * Nx)) + 1 #plus 1 because of 1-indexed arrays
#         lattice_points[idx] = i.*lattice.primitive_vectors[1,:] .+ j.*lattice.primitive_vectors[2,:]
#     end

#     return lattice_points
# end


function replicate_unit_cell(crystal::Crystal{D}, N::SVector{D}) where D
    @assert all(N .> 0) "Number of unit cells should be positive"

    #Probably a way to get LP an not allocate all this memory
    # lattice_points = get_lattice_points(crystal.lattice, N)
    N_atoms = prod(N) * length(crystal.basis)

    #Create flat arrays for atoms & coords
    atoms = SVector{N_atoms,Atom{D}}

    #Superimpose basis onto lattice points
    for i in range(N_atoms)
        n1 = 
        n2 =
        n3 =
        lattice_pt = crystal[n...]
        for basis_atom in crystal.basis         
            atoms[i] = Atom(basis_atom, lp .+ basis_atom.position)
        end
    end

    return atoms
end


#Monoatomic FCC

function FCC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), FaceCentered())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice, basis)
end



function BCC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), BodyCentered())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis)
end

# struct Diamond <: Crystal

# end

function Diamond(a, basis::Atom)

end