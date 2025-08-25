export
    Oblique,
    Rectangular, RectangularCentered,
    Square,
    Hexagonal


function Oblique(a::T, b::T, θ::T, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"q") where T
    lattice = BravaisLattice(ObliqueLattice(a, b, θ), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(b)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Oblique(a::T, b::T, θ::T, atomic_mass::Number, N::SVector{2}; 
                    charge = 0.0u"q", atomic_symbol::Symbol = :unknown) where T
    lattice = BravaisLattice(ObliqueLattice(a, b, θ), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(b)]; charge = charge, mass = atomic_mass)]
    return Crystal(lattice,basis,N)
end

#####

function Rectangular(a::T, b::T, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"q") where T
    lattice = BravaisLattice(RectangularLattice(a, b), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(b)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Rectangular(a::T, b::T, atomic_mass::Number, N::SVector{2}; 
                        charge = 0.0u"q", atomic_symbol::Symbol = :unknown) where T
    lattice = BravaisLattice(RectangularLattice(a, b), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(b)]; mass = atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#####

function RectangularCentered(a::T, b::T, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"q") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5) * a
    half_b = float_type(0.5) * b
    lattice = BravaisLattice(RectangularLattice(a, b), Primitive())
    basis = [Atom(atomic_symbol, [zero(a), zero(b)], charge = charge),
             Atom(atomic_symbol, [half_a, half_b], charge = charge)]
    return Crystal(lattice,basis,N)
end

function RectangularCentered(a::T, b::T, atomic_mass::Number, N::SVector{2}; 
                                charge = 0.0u"q", atomic_symbol::Symbol = :unknown) where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5) * a
    half_b = float_type(0.5) * b
    lattice = BravaisLattice(RectangularLattice(a, b), Primitive())
    basis = [Atom(atomic_symbol, [zero(a), zero(b)]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [half_a, half_b]; mass = atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#####

function Square(a::T, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"q") where T
    lattice = BravaisLattice(SquareLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Square(a::T, atomic_mass::Number, N::SVector{2}; 
                charge = 0.0u"q", atomic_symbol::Symbol = :unknown)  where T
    lattice = BravaisLattice(SquareLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a)]; mass = atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#####

function Hexagonal(a::T, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"q")  where T
    lattice = BravaisLattice(Hexagonal2DLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Hexagonal(a::T, atomic_mass::Number, N::SVector{2}; 
                    charge = 0.0u"q", atomic_symbol::Symbol = :unknown)  where T
    lattice = BravaisLattice(Hexagonal2DLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a)]; mass = atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end