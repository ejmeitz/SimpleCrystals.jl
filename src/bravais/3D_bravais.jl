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
function SC(a, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a), zero(a), zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function SC(a, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom([zero(a), zero(a), zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Monoatomic FCC -- 4 Atom Basis with SC Lattice Points
function FCC(a, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    T = typeof(ustrip(a))
    half_a = T(0.5)*a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, half_a, zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, zero(a), half_a], charge = charge),
             Atom(atomic_symbol, [zero(a), half_a, half_a], charge = charge)]
    return Crystal(lattice, basis, N)
end

function FCC(a, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    T = typeof(ustrip(a))
    half_a = T(0.5)*a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([half_a, half_a, zero(a)], atomic_mass, charge = charge),
             Atom([half_a, zero(a), half_a], atomic_mass, charge = charge),
             Atom([zero(a), half_a, half_a], atomic_mass, charge = charge)]
    return Crystal(lattice, basis, N)
end

# Monoatomic BCC -- 2 Atom Basis with SC Lattice Points
function BCC(a, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    T = typeof(ustrip(a))
    half_a = T(0.5)*a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [half_a, half_a, half_a], charge = charge)]
    return Crystal(lattice,basis,N)
end

function BCC(a, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    T = typeof(ustrip(a))
    half_a = T(0.5)*a
    lattice = BravaisLattice(CubicLattice(a), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([half_a, half_a, half_a], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#################
### Triclinic ###
#################

function Triclinic(a::T, b::T, c::T, α::T, β::T, γ::T, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(TriclinicLattice(a,b,c,α,β,γ), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Triclinic(a, b, c, α, β, γ, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(TriclinicLattice(a,b,c,α,β,γ), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

##################
### Monoclinic ###
##################

# Monoclinic Primitive
function Monoclinic(a, b, c, β, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(MonoclinicLattice(a,b,c,β), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Monoclinic(a, b, c, β, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(MonoclinicLattice(a,b,c,β), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Monoclinic Primitive w/ 2 Atom Basis
function MonoclinicBaseCentered(a, b, c, β, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(MonoclinicLattice(a,b,c,β), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [0.5*a, 0.5*b, zero(c)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function MonoclinicBaseCentered(a, b, c, β, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(MonoclinicLattice(a,b,c,β), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([0.5*a, 0.5*b, zero(c)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end
####################
### Orthorhombic ###
####################

function Orthorhombic(a, b, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(b),zero(c)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Orthorhombic(a, b, c, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom([zero(a),zero(b),zero(c)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Ortho Primitive w/ 2 Atom Basis
function OrthorhombicBaseCentered(a, b, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [0.5*a, 0.5*b, zero(c)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function OrthorhombicBaseCentered(a, b, c, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([0.5*a, 0.5*b, zero(c)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Ortho Primitive w/ 2 Atom Basis
function OrthorhombicBodyCentered(a, b, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [0.5*a, 0.5*b, 0.5*c], charge = charge)]
    return Crystal(lattice,basis,N)
end
function OrthorhombicBodyCentered(a, b, c, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([0.5*a, 0.5*b, 0.5*c], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Ortho Primitive w/ 4 Atom Basis
function OrthorhombicFaceCentered(a, b, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [0.5*a, 0.5*b, zero(c)], charge = charge),
             Atom(atomic_symbol, [0.5*a, zero(b), 0.5*c], charge = charge),
             Atom(atomic_symbol, [zero(a), 0.5*b, 0.5*c], charge = charge)]
    return Crystal(lattice,basis,N)
end

function OrthorhombicFaceCentered(a, b, c, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(OrthorhombicLattice(a,b,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([0.5*a, 0.5*b, zero(c)], atomic_mass, charge = charge),
             Atom([0.5*a, zero(b), 0.5*c], atomic_mass, charge = charge),
             Atom([zero(a), 0.5*b, 0.5*c], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

##################
### Tetragonal ###
##################

# Tetragonal Primitive
function Tetragonal(a, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(TetragonalLattice(a,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(c)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Tetragonal(a, c, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(TetragonalLattice(a,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(c)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

# Tetragonal Primitive w/ 2 Atom Basis
function TetragonalBodyCentered(a, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(TetragonalLattice(a,c), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge),
             Atom(atomic_symbol, [0.5*a, 0.5*a, 0.5*c], charge = charge)]
    return Crystal(lattice,basis,N)
end

function TetragonalBodyCentered(a, c, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(TetragonalLattice(a,c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge),
             Atom([0.5*a, 0.5*a, 0.5*c], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

#################
### Hexagonal ###
#################

function Rhombohedral(a, α, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(RhombohedralLattice(a, α), Primitive())
    basis = [Atom(atomic_symbol, [zero(a),zero(a),zero(a)], charge = charge)]
    return Crystal(lattice,basis,N)
end

function Rhombohedral(a, α, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(RhombohedralLattice(a, α), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end

function Hexagonal(a, c, atomic_mass::Number, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(HexagonalLattice(a, c), Primitive())
    basis = [Atom([zero(a),zero(a),zero(a)], atomic_mass, charge = charge)]
    return Crystal(lattice,basis,N)
end