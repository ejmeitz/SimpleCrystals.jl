export
    Atom

#Load periodic table data
periodic_table = PeriodicTable.elements

struct Atom{D,C,M}
    sym::Symbol
    position::SVector{D}
    charge::C
    mass::M
end

function Atom(sym, position; charge =0.0u"C", mass = periodic_table[sym].atomic_mass)
    return Atom{length(position),typeof(charge),typeof(mass)}(sym, position, charge, mass)
end

#Construct atom from another, but with new position
function Atom(atom::Atom{D,C,M}, position::SVector{D}) where {D,C,M}
    return Atom{D,C,M}(atom.sym, position, atom.charge, atom.mass)
end

AtomsBase.atomic_symbol(atom::Atom) = atom.sym
AtomsBase.atomic_mass(atom::Atom) = atom.mass
AtomsBase.atomic_number(atom::Atom) = periodic_table[atom.sym].number
AtomsBase.position(atom::Atom) = atom.position