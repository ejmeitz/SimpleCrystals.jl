export Honeycomb

# Hexagonal 2D Lattice -- 2 atom basis
function Honeycomb(a, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"C")
    d = a*sqrt(3)
    lattice = BravaisLattice(Hexagonal2DLattice(a*sqrt(3)), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [0.5*d, 0.5*a], charge = charge),]
    return Crystal(lattice,basis,N)
end

function Honeycomb(a, atomic_mass::Number, N::SVector{2}; charge = 0.0u"C")
    d = a*sqrt(3)
    lattice = BravaisLattice(Hexagonal2DLattice(a*sqrt(3)), Primitive())
    basis = [Atom([zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([0.5*d, 0.5*a], atomic_mass, charge = charge),]
    return Crystal(lattice,basis,N)
end