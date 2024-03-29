export to_xyz

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