export to_xyz, to_ucposcar, to_ssposcar

function to_xyz(crystal::Crystal{3}, outpath)
    N_atoms = length(crystal)
    
    file = open(outpath, "w")
    println(file,string(N_atoms))
    println(file,"Comment Line")
    for atom in crystal
        println(file,"$(string(atomic_symbol(atom))) $(ustrip(atom.position[1])) $(ustrip(atom.position[2])) $(ustrip(atom.position[3]))")
    end
    close(file)
end

function to_xyz(crystal::Crystal{2}, outpath)
    N_atoms = length(crystal)
    
    file = open(outpath, "w")
    println(file,string(N_atoms))
    println(file,"Comment Line")
    for atom in crystal
        println(file,"$(string(atomic_symbol(atom))) $(ustrip(atom.position[1])) $(ustrip(atom.position[2]))")
    end
    close(file)
end

"""
    cartesian_to_fractional(latt::BravaisLattice{3}, x::AbstractVector) where T

Convert a Cartesian coordinate `x` into fractional coordinates w.r.t. the
primitive vectors of `latt`.
"""
cartesian_to_fractional(lattice_vectors::AbstractMatrix, x::AbstractVector) = lattice_vectors \ x


function convert_to_ang_and_strip(x)
    if eltype(x) <: Unitful.Length
        x = ustrip.(u"Å", x)
    else
        @error "Expected length units, got $(unit(first(x)))"
    end
end

function write_poscar(cryst::Crystal{D}, cell::AbstractMatrix, path::AbstractString) where D
    # 1) collect element order & counts
    all_syms = [atomic_symbol(atom) for atom in cryst.atoms]
    # preserve order of first appearance in basis
    order = unique([atomic_symbol(atom) for atom in cryst.basis])
    counts = [count(==(sym), all_syms) for sym in order]

    open(path, "w") do io
        println(io, "POSCAR")                   # comment line
        println(io, "1.0")                      # universal scale

        # primitive vectors as rows
        for i in 1:D
            v = convert_to_ang_and_strip(cell[:,i])
            for el in v
                @printf(io, "%.15f ", el)
            end
            println(io)
        end

        # elements and counts
        println(io, join(order, " "))
        println(io, join(counts, " "))

        println(io, "Direct")  # fractional coords

        # write atoms in block order
        for sym in order
            for atom in cryst.atoms
                if atomic_symbol(atom) == sym
                    x = convert_to_ang_and_strip(position(atom))
                    vec = convert_to_ang_and_strip(cell)
                    x_frac = cartesian_to_fractional(vec, x)
                    for el in x_frac
                        @printf(io, "%.15f ", el)
                    end
                    println(io)
                end
            end
        end
    end
end

"""
to_ucposcar(cryst::Crystal{3}; filename="ucposcar")

Automatically extracts the *unit cell* from `cryst` (i.e. N_unit_cells ≡ (1,1,1))
and writes it to a POSCAR-style file named `filename`.
"""
function to_ucposcar(cryst::Crystal{D}, filename::AbstractString="ucposcar") where D
    
    n_cells = @SVector ones(Int, D)

    uc = Crystal(
        cryst.lattice,
        cryst.basis,
        n_cells,
        cryst.basis
    )
    write_poscar(uc, cryst.lattice.primitive_vectors, filename)
end

"""
    supercell(cryst::Crystal{3}; filename="ssposcar")

Writes the full supercell defined by `cryst.N_unit_cells` to a
POSCAR-style file named `filename`.
"""
function to_ssposcar(cryst::Crystal, filename::AbstractString="ssposcar")
    cell = cryst.lattice.primitive_vectors .* cryst.N_unit_cells
    write_poscar(cryst, cell, filename)
end