export Honeycomb

# Hexagonal 2D Lattice -- 2 atom basis
function Honeycomb(a::T, atomic_symbol::Symbol, N::SVector{2}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    d = a*sqrt(float_type(3.0))
    half_a = float_type(0.5) * a
    half_d = float_type(0.5) * d
    lattice = BravaisLattice(Hexagonal2DLattice(a*sqrt(3)), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_d, half_a], charge = charge),]
    return Crystal(lattice,basis,N)
end

function Honeycomb(a::T, atomic_mass::Number, N::SVector{2}; 
                    charge = 0.0u"C", atomic_symbol::Symbol = :unknown) where T
    float_type = typeof(ustrip(a))
    d = a*sqrt(float_type(3.0))
    half_a = float_type(0.5) * a
    half_d = float_type(0.5) * d
    lattice = BravaisLattice(Hexagonal2DLattice(a*sqrt(3)), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a)]; mass = atomic_mass, charge = charge),
             Atom(atomic_symbol, [half_d, half_a], mass = atomic_mass, charge = charge),]
    return Crystal(lattice,basis,N)
end