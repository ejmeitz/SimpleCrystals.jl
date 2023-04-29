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
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a), zero(a), zero(a)), charge = charge)]
    return Crystal(lattice,basis,N)
end

# Monoatomic FCC -- 4 Atom Basis with SC Lattice Points
function FCC(a, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(0.5*a, 0.5*a, zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(zero(a), 0.5*a, 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, zero(a), 0.5*a), charge = charge),
             Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice, basis,N)
end

# Monoatomic BCC -- 2 Atom Basis with SC Lattice Points
function BCC(a, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Cubic(a), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*a, 0.5*a), charge = charge)]
    return Crystal(lattice,basis,N)
end

#################
### Triclinic ###
#################

function Triclinic(a, b, c, α, β, γ, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Triclinic(a,b,c,α,β,γ), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis,N)
end

##################
### Monoclinic ###
##################

function Monoclinic(a, b, c, β, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Monoclinic(a,b,c,β), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis,N)
end

# Monoclinic Primitive w/ 2 Atom Basis
function MonoclinicBaseCentered(a, b, c, β, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Monoclinic(a,b,c,β), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*b, zero(c)), charge = charge)]
    return Crystal(lattice,basis,N)
end

####################
### Orthorhombic ###
####################

function Orthorhombic(a, b, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Orthorhombic(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(b),zero(c)), charge = charge)]
    return Crystal(lattice,basis,N)
end

# Ortho Primitive w/ 2 Atom Basis
function OrthorhombicBaseCentered(a, b, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Orthorhombic(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*b, zero(c)), charge = charge)]
    return Crystal(lattice,basis,N)
end

# Ortho Primitive w/ 2 Atom Basis
function OrthorhombicBodyCentered(a, b, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Orthorhombic(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*b, 0.5*c), charge = charge)]
    return Crystal(lattice,basis,N)
end

# Ortho Primitive w/ 4 Atom Basis
function OrthorhombicFaceCentered(a, b, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Orthorhombic(a,b,c), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*b, zero(c)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, zero(b), 0.5*c), charge = charge),
             Atom(atomic_symbol, SVector(zero(a), 0.5*b, 0.5*c), charge = charge)]
    return Crystal(lattice,basis,N)
end

##################
### Tetragonal ###
##################

# Ortho Primitive w/ 4 Atom Basis
function Tetragonal(a, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Tetragonal(a,c), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(c)), charge = charge)]
    return Crystal(lattice,basis,N)
end

# Tetragonal Primitive w/ 2 Atom Basis
function TetragonalBodyCentered(a, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Tetragonal(a,c), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge),
             Atom(atomic_symbol, SVector(0.5*a, 0.5*a, 0.5*c), charge = charge)]
    return Crystal(lattice,basis,N)
end

#################
### Hexagonal ###
#################

function Rhombohedral(a, α, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Rhombohedral(a, α), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis,N)
end

function Hexagonal(a, c, atomic_symbol::Symbol, N::SVector{3}; charge = 0.0u"C")
    lattice = BravaisLattice(Hexagonal(a, c), Primitive())
    basis = [Atom(atomic_symbol, SVector(zero(a),zero(a),zero(a)), charge = charge)]
    return Crystal(lattice,basis,N)
end