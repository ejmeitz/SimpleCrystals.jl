export
    Oblique,
    Rectangular, RectangularCentered,
    Square,
    Hexagonal


function Oblique(a, b, θ, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Oblique(a, b, θ), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(b)), charge = charge)]
    return Crystal(lattice,basis)
end

function Rectangular(a, b, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Rectangular(a, b), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(b)), charge = charge)]
    return Crystal(lattice,basis)
end

function RectangularCentered(a, b, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Rectangular(a, b), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(b)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*b), charge = charge)]
    return Crystal(lattice,basis)
end

function Square(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Square(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis)
end

function Hexagonal(a, atomic_symbol::Symbol; charge = 0.0u"C")
    lattice = BravaisLattice(Hexagonal2D(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis)
end