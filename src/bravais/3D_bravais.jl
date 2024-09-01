export
    SC, FCC, BCC,
    Triclinic,
    Monoclinic, MonoclinicBaseCentered,
    Orthorhombic,OrthorhombicBaseCentered, OrthorhombicBodyCentered, OrthorhombicFaceCentered,
    Tetragonal, TetragonalBodyCentered,
    Rhombohedral, Hexagonal

#############
### Cubic ###
#############

# Monoatomic SC -- 1 Atom Basis with SC Lattice Points
function SC(a::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a), zero(a), zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function SC(a::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom([zero(a), zero(a), zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Monoatomic FCC -- 4 Atom Basis with SC Lattice Points
function FCC(a::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, half_a, zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, zero(a), half_a], charge = charge),
             Atom(atomic_symbol, [zero(a), half_a, half_a], charge = charge)]
    return Crystal(lattice, basis, N)
end

function FCC(a, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([half_a, half_a, zero(a)], atomic_mass, charge = charge),
             Atom([half_a, zero(a), half_a], atomic_mass, charge = charge),
             Atom([zero(a), half_a, half_a], atomic_mass, charge = charge)]
    return Crystal(lattice, basis, N)
end

# Monoatomic BCC -- 2 Atom Basis with SC Lattice Points
function BCC(a::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, half_a, half_a], charge = charge)]
    return Crystal(lattice,basis,N)
end

function BCC(a::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([half_a, half_a, half_a], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#################
### Triclinic ###
#################

function Triclinic(a::T, b::T, c::T, α::T, β::T, γ::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(TriclinicLattice(a,b,c,α,β,γ), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Triclinic(a::T, b::T, c::T, α::T, β::T, γ::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(TriclinicLattice(a,b,c,α,β,γ), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

##################
### Monoclinic ###
##################

# Monoclinic Primitive
function Monoclinic(a::T, b::T, c::T, β::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(MonoclinicLattice(a,b,c,β), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Monoclinic(a::T, b::T, c::T, β::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(MonoclinicLattice(a,b,c,β), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Monoclinic Primitive w/ 2 Atom Basis
function MonoclinicBaseCentered(a::T, b::T, c::T, β::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_b = float_type(0.5)*b
    lattice = BravaisLattice(MonoclinicLattice(a,b,c,β), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, half_b, zero(c)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function MonoclinicBaseCentered(a::T, b::T, c::T, β::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_b = float_type(0.5)*b
    lattice = BravaisLattice(MonoclinicLattice(a,b,c,β), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([half_a, half_b, zero(c)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end
####################
### Orthorhombic ###
####################

function Orthorhombic(a::T, b::T, c::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(b),zero(c)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Orthorhombic(a::T, b::T, c::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom([zero(a),zero(b),zero(c)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Ortho Primitive w/ 2 Atom Basis
function OrthorhombicBaseCentered(a::T, b::T, c::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_b = float_type(0.5)*b
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, half_b, zero(c)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function OrthorhombicBaseCentered(a::T, b::T, c::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_b = float_type(0.5)*b
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([half_a, half_b, zero(c)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Ortho Primitive w/ 2 Atom Basis
function OrthorhombicBodyCentered(a::T, b::T, c::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_b = float_type(0.5)*b
    half_c = float_type(0.5)*c
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, half_b, half_c], charge = charge)]
    return Crystal(lattice,basis,N)
end
function OrthorhombicBodyCentered(a::T, b::T, c::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_b = float_type(0.5)*b
    half_c = float_type(0.5)*c
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([half_a, half_b, half_c], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Ortho Primitive w/ 4 Atom Basis
function OrthorhombicFaceCentered(a::T, b::T, c::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_b = float_type(0.5)*b
    half_c = float_type(0.5)*c
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, half_b, zero(c)], charge = charge),
             Atom(atomic_symbol, [half_a, zero(b), half_c], charge = charge),
             Atom(atomic_symbol, [zero(a), half_b, half_c], charge = charge)]
    return Crystal(lattice,basis,N)
end

function OrthorhombicFaceCentered(a::T, b::T, c::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_b = float_type(0.5)*b
    half_c = float_type(0.5)*c
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([half_a, half_b, zero(c)], atomic_mass, charge = charge),
             Atom([half_a, zero(b), half_c], atomic_mass, charge = charge),
             Atom([zero(a), half_b, half_c], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

##################
### Tetragonal ###
##################

# Tetragonal Primitive
function Tetragonal(a::T, c::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(TetragonalLattice(a,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(c)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Tetragonal(a::T, c::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(TetragonalLattice(a,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(c)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Tetragonal Primitive w/ 2 Atom Basis
function TetragonalBodyCentered(a::T, c::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_c = float_type(0.5)*c
    lattice = BravaisLattice(TetragonalLattice(a,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, half_a, half_c], charge = charge)]
    return Crystal(lattice,basis,N)
end

function TetragonalBodyCentered(a::T, c::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    float_type = typeof(ustrip(a))
    half_a = float_type(0.5)*a
    half_c = float_type(0.5)*c
    lattice = BravaisLattice(TetragonalLattice(a,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([half_a, half_a, half_c], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#################
### Hexagonal ###
#################

function Rhombohedral(a::T, α::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(RhombohedralLattice(a, α), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Rhombohedral(a::T, α::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(RhombohedralLattice(a, α), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

function Hexagonal(a::T, c::T, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C") where T
    lattice = BravaisLattice(HexagonalLattice(a, c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end