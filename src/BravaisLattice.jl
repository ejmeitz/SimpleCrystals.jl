export
    Crystal, BravaisLattice,
    CubicLattice, SquareLattice, OrthorhombicLattice, RectangularLattice,
    MonoclinicLattice, ObliqueLattice, TriclinicLattice, TetragonalLattice,
    RhombohedralLattice, HexagonalLattice, Hexagonal2DLattice,
    Primitive, FaceCentered, BodyCentered, BaseCentered, Centered
    #some of the exported BravaisLattice's have the same name as the exported crystals, but they take different parameters
 
#########
# TYEPS #
#########

abstract type CrystalFamily{D} end
abstract type CenteringType end

#Centering types for multiple dispatch
struct Primitive <: CenteringType end
struct FaceCentered <: CenteringType end
struct BodyCentered <: CenteringType end
struct BaseCentered <: CenteringType end
struct Centered <: CenteringType end #2D


struct BravaisLattice{D} #D is dimension
    crystal_family::CrystalFamily{D}
    centering_type::CenteringType
    primitive_vectors::MMatrix{D,D}
end

function BravaisLattice(cf::CrystalFamily{D}, ct::CenteringType) where D
    p_vec = get_primitive_vectors(cf,ct)
    return BravaisLattice{D}(cf, ct, p_vec)
end

#Returns R = n1*a1 + n2*a2 + n3*a3, where a are the primitive lattice vectors
function Base.getindex(lattice::BravaisLattice{D}, indices::Vararg{Integer,D}) where D
    return SVector{D}(sum(indices .* lattice.primitive_vectors, dims = 1))
end

#####################################################

struct CubicLattice{LC} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
end

CubicLattice(a) = CubicLattice{typeof(a)}(SVector(a,a,a))

struct SquareLattice{LC} <: CrystalFamily{2}
    lattice_constants::SVector{2,LC}
end

SquareLattice(a) = SquareLattice{typeof(a)}(SVector(a,a))

#####################################################

struct OrthorhombicLattice{LC} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
end

OrthorhombicLattice(a,b,c) = OrthorhombicLattice{typeof(a)}(SVector(a,b,c))


struct RectangularLattice{LC} <: CrystalFamily{2}
    lattice_constants::SVector{2,LC}
end

RectangularLattice(a,b) = RectangularLattice{typeof(a)}(SVector(a,b))

#####################################################

struct MonoclinicLattice{LC,LA} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

MonoclinicLattice(a, b, c, β) = MonoclinicLattice{typeof(a),typeof(β)}(SVector(a,b,c), SVector(90u"°", β, 90u"°"))

# 2D equivalent
struct ObliqueLattice{LC,LA} <: CrystalFamily{2}
    lattice_constants::SVector{2,LC}
    lattice_angle::LA
end

ObliqueLattice(a, b, θ) = ObliqueLattice{typeof(a),typeof(θ)}(SVector(a,b), θ)
#####################################################

struct TriclinicLattice{LC,LA} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

TriclinicLattice(a, b, c, α, β, γ) = TriclinicLattice{typeof(a),typeof(β)}(SVector(a,b,c), SVector(α, β, γ))

#####################################################

struct TetragonalLattice{LC} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
end

TetragonalLattice(a, c) = TetragonalLattice{typeof(a)}(SVector(a,a,c))

#####################################################

struct RhombohedralLattice{LC,LA} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

function RhombohedralLattice(a, α)
    return RhombohedralLattice{typeof(a),typeof(α)}(SVector(a,a,a),SVector(α, α, α))
end

#####################################################

struct HexagonalLattice{LC,LA} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

#Defined with 60deg instead of 120 so that domain is triclinic
HexagonalLattice(a,c) = HexagonalLattice{typeof(a),typeof(90u"°")}(SVector(a,a,c),SVector(90u"°", 90u"°", 60u"°"))

struct Hexagonal2DLattice{LC,LA} <: CrystalFamily{2}
    lattice_constants::SVector{2,LC}
    lattice_angle::LA
end

Hexagonal2DLattice(a) = Hexagonal2DLattice{typeof(a),typeof(120u"°")}(SVector(a,a), 120u"°")

#####################################################

function get_primitive_vectors(cf::CrystalFamily{3}, ct::Primitive)
    primitive_vectors = MMatrix{3,3}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)

    if hasfield(typeof(cf), :lattice_angles)
        rotate_primitive_vectors!(primitive_vectors, cf.lattice_constants[2], cf.lattice_constants[3],
             cf.lattice_angles[1], cf.lattice_angles[2], cf.lattice_angles[3])
    end

    return primitive_vectors
end

FaceCenteredSupportedTypes = Union{CubicLattice, OrthorhombicLattice}
function get_primitive_vectors(cf::FaceCenteredSupportedTypes, ct::FaceCentered)
    primitive_vectors = MMatrix{3,3}([0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]) 
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)

    if hasfield(typeof(cf), :lattice_angles)
        rotate_primitive_vectors!(primitive_vectors, cf.lattice_constants[2], cf.lattice_constants[3],
             cf.lattice_angles[1], cf.lattice_angles[2], cf.lattice_angles[3])
    end

    return primitive_vectors
end

BodyCenteredSupportedTypes = Union{CubicLattice, OrthorhombicLattice, TetragonalLattice}
function get_primitive_vectors(cf::BodyCenteredSupportedTypes, ct::BodyCentered)
    primitive_vectors = MMatrix{3,3}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.5  0.5  0.5])
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)

    if hasfield(typeof(cf), :lattice_angles)
        rotate_primitive_vectors!(primitive_vectors, cf.lattice_constants[2], cf.lattice_constants[3],
             cf.lattice_angles[1], cf.lattice_angles[2], cf.lattice_angles[3])
    end

    return primitive_vectors
end

BaseCenteredSupportedTypes = Union{MonoclinicLattice, OrthorhombicLattice}
function get_primitive_vectors(cf::BaseCenteredSupportedTypes, ct::BaseCentered)
    primitive_vectors = MMatrix{3,3}([1.0 1.0 0.0; 1.0 -1.0 0.0; 0.0  0.0  1.0])
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)

    if hasfield(typeof(cf), :lattice_angles)
        rotate_primitive_vectors!(primitive_vectors, cf.lattice_constants[2], cf.lattice_constants[3],
             cf.lattice_angles[1], cf.lattice_angles[2], cf.lattice_angles[3])
    end

    return primitive_vectors
end

function get_primitive_vectors(cf::CrystalFamily{2}, ct::Primitive)
    primitive_vectors = MMatrix{2,2}([1.0 0.0; 0.0 1.0])
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)
    
    if hasfield(typeof(cf), :lattice_angle)
        β = cf.lattice_angle - 90u"°"
        primitive_vectors[2,:] = [cos(β) -sin(β); sin(β)  cos(β)] * primitive_vectors[2,:]
    end
    return primitive_vectors
end

function get_primitive_vectors(cf::RectangularLattice, ct::Centered)
    primitive_vectors = MMatrix{2,2}([0.5 0.5; 0.5 -0.5])
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)
    return primitive_vectors
end


#####################################################
# Helper functions

# For arbitrary rotaions
# Maintains right-handed coordinate system
function rotate_primitive_vectors!(p_vec, b, c, α, β, γ)
    # Rotate vector along y-axis
    p_vec[2,:] = [b*cos(γ), b*sin(γ), zero(b)]
    # Rotate vector along z-axis
    cx = c*cos(β)
    cy = c*(cos(α) - (cos(β)*cos(γ)))/sin(γ)
    cz = sqrt(c^2 - cx^2 - cy^2)
    p_vec[3,:] = [cx, cy, cz]
    return p_vec
end
