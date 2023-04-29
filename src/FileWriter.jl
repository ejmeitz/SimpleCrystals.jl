export to_xyz

function to_xyz(crystal::Crystal{3}, outpath)
    atoms = get_coordinates(crystal)
    
    file = open(outpath, "w")
    println(file,string(length(atoms)))
    println(file,"Comment Line")
    for atom in atoms
        println(file,"$(string(atom.sym)) $(ustrip(atom.position[1])) $(ustrip(atom.position[2])) $(ustrip(atom.position[3]))")
    end
    close(file)
end

function to_xyz(crystal::Crystal{2}, outpath)
    atoms = get_coordinates(crystal)
    
    file = open(outpath, "w")
    println(file,string(length(atoms)))
    println(file,"Comment Line")
    for atom in atoms
        println(file,"$(string(atom.sym)) $(ustrip(atom.position[1])) $(ustrip(atom.position[2]))")
    end
    close(file)
end