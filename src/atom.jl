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

function Atom(sym::Symbol, position; charge =0.0u"C", mass = periodic_table[sym].atomic_mass)
    return Atom{length(position),typeof(charge),typeof(mass)}(sym, position, charge, mass)
end

AtomsBase.atomic_symbol(atom::Atom) = atom.sym
AtomsBase.atomic_mass(atom::Atom) = atom.mass
AtomsBase.atomic_number(atom::Atom) = periodic_table[atom.sym].number
AtomsBase.position(atom::Atom) = atom.position
AtomsBase.n_dimensions(atom::Atom) = length(atom.position)