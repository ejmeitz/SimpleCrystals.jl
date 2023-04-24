export
    SC, FCC, BCC,
    Rhombohedral

# These strucrtures use the cubic conventional cell so that a
# cubic simulation cell can be used.

# Monoatomic SC -- 1 Atom Basis with SC Lattice Points
function SC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a), zero(a), zero(a)), charge = charge)]
    return Crystal(lattice,basis)
end

# Monoatomic FCC -- 4 Atom Basis with SC Lattice Points
function FCC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(0.5*a, 0.5*a, zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(zero(a), 0.5*a, 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, zero(a), 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice, basis)
end

# Monoatomic BCC -- 2 Atom Basis with SC Lattice Points
function BCC(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*a, 0.5*a), charge = charge)]
    return Crystal(lattice,basis)
end

# All the strucrtures below require a triclinic box when imported to a MD code

# Rhombohedral lattice -- 1 atom basis
# Simple cubic with α not equal to 90
function Rhombohedral(a, α, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Rhombohedral(a, α), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis)
end