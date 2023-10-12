export
    Atom

#Load periodic table data
periodic_table = PeriodicTable.elements

struct Atom{D,T,C,M}
    atomic_symbol::Symbol
    position::Vector{T}
    charge::C
    atomic_mass::M
end

function Atom(sym::Symbol, position; charge =0.0u"q", mass = periodic_table[sym].atomic_mass)
    return Atom{length(position),eltype(position),typeof(charge),typeof(mass)}(sym, position, charge, mass)
end

function Atom(position, mass::Number; charge = 0.0u"q")
    return Atom{length(position),eltype(position),typeof(charge), typeof(mass)}(:unknown, position, charge, mass)
end


charge(atom::Atom) = atom.charge

Base.keys(atom::Atom) = fieldnames(typeof(atom))
Base.haskey(atom::Atom, x::Symbol) = hasfield(Atom, x)
Base.getindex(atom::Atom, x::Symbol) = hasfield(Atom, x) ? getfield(atom, x) : error("No field `$x` in Atom object. Allowed keys are $(keys(atom)).")
Base.get(atom::Atom, x::Symbol, default) = hasfield(Atom, x) ? getfield(atom,x) : default
Base.pairs(atom::Atom) = (k => atom[k] for k in keys(atom))

AtomsBase.atomic_symbol(atom::Atom) = atom.atomic_symbol
AtomsBase.atomic_mass(atom::Atom) = atom.atomic_mass
AtomsBase.atomic_number(atom::Atom) = (atom.atomic_symbol == :unknown) ? :unknown : periodic_table[atom.atomic_symbol].number
AtomsBase.position(atom::Atom) = atom.position
AtomsBase.velocity(atom::Atom) = missing
AtomsBase.n_dimensions(::Atom{D}) where D = D

function Base.show(io::IO, atom::Atom)
    print(io, "Atom at $(round.(ustrip.(atom.position), digits = 3)), with charge: $(charge(atom)) and mass : $(atomic_mass(atom)) ")
end


function Base.show(io::IO, ::MIME"text/plain", atom::Atom)
    print(io, "Atom at $(round.(ustrip.(atom.position), digits = 3)), with charge: $(charge(atom)) and mass : $(atomic_mass(atom)) ")
end
