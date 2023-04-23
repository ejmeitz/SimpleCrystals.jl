function to_xyz(crystal::Crystal{D}, N::SVector{D}, outpath)
    atoms = replicate_unit_cell(crystal::Crystal{D}, N::SVector{D})
    
    file = open(outpath, "w")
    println(file,string(length(atoms)))
    println(file,"Comment Line")
    for atom in atoms
        println(file,"$(ustrip(atom.position[1])) $(ustrip(atom.position[2])) $(ustrip(atom.position[3]))")
    end
    close(file)
end