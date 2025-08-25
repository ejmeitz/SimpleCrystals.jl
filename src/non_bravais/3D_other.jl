export
    Diamond,
    HCP
    
# This implementation of diamond uses a cubic conventional cell so that a
# cubic simulation cell can be used.

# Monoatomic Diamond -- 8 Atom Basis with SC Lattice Points
function Diamond(a::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"q") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5) * a
    quarter_a = float_type(0.25) * a
    three_quarter_a = float_type(0.75) * a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [zero(a), half_a, half_a], charge = charge),
             Atom(atomic_symbol, [half_a, zero(a), half_a], charge = charge),
             Atom(atomic_symbol, [half_a, half_a, zero(a)], charge = charge),
             Atom(atomic_symbol, [quarter_a, quarter_a, quarter_a], charge = charge),
             Atom(atomic_symbol, [quarter_a, three_quarter_a, three_quarter_a], charge = charge),
             Atom(atomic_symbol, [three_quarter_a, quarter_a, three_quarter_a], charge = charge),
             Atom(atomic_symbol, [three_quarter_a, three_quarter_a, quarter_a], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Diamond(a::T, atomic_mass::Number, N::SVector{3}; 
                    charge = 0.0u"q", atomic_symbol::Symbol = :unknown) where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5) * a
    quarter_a = float_type(0.25) * a
    three_quarter_a = float_type(0.75) * a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [zero(a), half_a, half_a]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [half_a, zero(a), half_a]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [half_a, half_a, zero(a)]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [quarter_a, quarter_a, quarter_a]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [quarter_a, three_quarter_a, three_quarter_a]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [three_quarter_a, quarter_a, three_quarter_a]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [three_quarter_a, three_quarter_a, quarter_a]; mass = atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# HCP requires a triclinic simulation box when imported to a MD code

# Hexagonal Lattice with a = c -- 2 atom basis
function HCP(a::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"q") where T
    float_type = typeof(ustrip(a))
    c = a*sqrt(float_type(8/3))
    half_a = float_type(0.5) * a
    half_c = float_type(0.5) * c
    b = float_type(a/(2*sqrt(3)))
    lattice = BravaisLattice(HexagonalLattice(a,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, b, half_c], charge = charge)]
    return Crystal(lattice,basis,N)
end

function HCP(a::T, atomic_mass::Number, N::SVector{3}; 
                charge = 0.0u"q", atomic_symbol::Symbol = :unknown) where T
    float_type = typeof(ustrip(a))
    c = a*sqrt(float_type(8/3))
    half_a = float_type(0.5) * a
    half_c = float_type(0.5) * c
    b = float_type(a/(2*sqrt(3)))
    lattice = BravaisLattice(HexagonalLattice(a,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [half_a, b, half_c]; mass = atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end