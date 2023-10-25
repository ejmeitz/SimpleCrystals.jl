export
    Oblique,
    Rectangular, RectangularCentered,
    Square,
    Hexagonal


function Oblique(a, b, θ, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(ObliqueLattice(a, b, θ), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(b)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Oblique(a, b, θ, atomic_mass::Number, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(ObliqueLattice(a, b, θ), Primitive())
    basis = [Atom([zero(a),zero(b)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#####

function Rectangular(a, b, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(RectangularLattice(a, b), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(b)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Rectangular(a, b, atomic_mass::Number, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(RectangularLattice(a, b), Primitive())
    basis = [Atom([zero(a),zero(b)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#####

function RectangularCentered(a, b, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(RectangularLattice(a, b), Primitive())
    basis = [Atom(atomic_symbol, [zero(a), zero(b)], charge = charge),
             Atom(atomic_symbol, [0.5*a, 0.5*b], charge = charge)]
    return Crystal(lattice,basis,N)
end

function RectangularCentered(a, b, atomic_mass::Number, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(RectangularLattice(a, b), Primitive())
    basis = [Atom([zero(a), zero(b)], atomic_mass, charge = charge),
             Atom([0.5*a, 0.5*b], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#####

function Square(a, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(SquareLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Square(a, atomic_mass::Number, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(SquareLattice(a), Primitive())
    basis = [Atom([zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#####

function Hexagonal(a, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(Hexagonal2DLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Hexagonal(a, atomic_mass::Number, N::SVector{2}; charge = 0.0u"C")
    lattice = BravaisLattice(Hexagonal2DLattice(a), Primitive())
    basis = [Atom([zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end