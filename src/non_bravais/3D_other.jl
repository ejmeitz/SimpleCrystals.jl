export
    Diamond,
    HCP
    
# This implementation of diamond uses a cubic conventional cell so that a
# cubic simulation cell can be used.

# Monoatomic Diamond -- 8 Atom Basis with SC Lattice Points
function Diamond(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(zero(a), 0.5*a, 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, zero(a), 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*a, zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.25*a, 0.25*a, 0.25*a), charge = charge),
             Atom(atomic_symbol, SVector(0.25*a, 0.75*a, 0.75*a), charge = charge),
             Atom(atomic_symbol, SVector(0.75*a, 0.25*a, 0.75*a), charge = charge),
             Atom(atomic_symbol, SVector(0.75*a, 0.75*a, 0.25*a), charge = charge)]
    return Crystal(lattice,basis)
end

# HCP requires a triclinic simulation box when imported to a MD code

# Hexagonal Lattice with a = c -- 2 atom basis
function HCP(a, atomic_symbol::Symbol; charge = 0.0u"C")
    c = a*sqrt(8/3)
    lattice = BravaisLattice(Hexagonal(a,c), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, a/(2*sqrt(3)), 0.5*c), charge = charge)]
    return Crystal(lattice,basis)
end