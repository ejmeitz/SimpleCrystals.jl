export
    Crystal
 
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


#####################################################

struct Cubic{LC} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
end

Cubic(a) = Cubic{typeof(a)}(SVector(a,a,a))

struct Square{LC} <: CrystalFamily{2}
    lattice_constants::SVector{2,LC}
end

Square(a) = Square{typeof(a)}(SVector(a,a))

#####################################################

struct Orthorhombic{LC} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
end

Orthorhombic(a,b,c) = Orthorhombic{typeof(a)}(SVector(a,b,c))


struct Rectangular{LC} <: CrystalFamily{2}
    lattice_constants::SVector{2,LC}
end

Rectangular(a,b) = Rectangular{typeof(a)}(SVector(a,b))

#####################################################

struct Monoclinic{LC,LA} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

Monoclinic(a, b, c, β) = Monoclinic{typeof(a),typeof(β)}(SVector(a,b,c), SVector(90u"°", β, 90u"°"))

# 2D equivalent
struct Oblique{LC,LA} <: CrystalFamily{2}
    lattice_constants::SVector{2,LC}
    lattice_angle::LA
end

Oblique(a, b, θ) = Oblique{typeof(a),typeof(θ)}(SVector(a,b), θ)
#####################################################

struct Triclinic{LC,LA} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

Triclinic(a, b, c, α, β, γ) = Triclinic{typeof(a),typeof(β)}(SVector(a,b,c), SVector(α, β, γ))

#####################################################

struct Tetragonal{LC} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
end

Tetragonal(a, c) = Tetragonal{typeof(a)}(SVector(a,a,c))

#####################################################

struct Rhombohedral{LC,LA} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

function Rhombohedral(a, α)
    return Rhombohedral{typeof(a),typeof(α)}(SVector(a,a,a),SVector(α, α, α))
end

#####################################################

struct Hexagonal{LC,LA} <: CrystalFamily{3}
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

Hexagonal(a,c) = Hexagonal{typeof(a),typeof(90u"°")}(SVector(a,a,c),SVector(90u"°", 90u"°", 120u"°"))

struct Hexagonal2D{LC,LA} <: CrystalFamily{2}
    lattice_constants::SVector{2,LC}
    lattice_angle::LA
end

Hexagonal2D(a) = Hexagonal2D{typeof(a),typeof(120u"°")}(SVector(a,a), 120u"°")

#####################################################

function get_primitive_vectors(cf::CrystalFamily{3}, ct::Primitive)
    primitive_vectors = MMatrix{3,3}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)

    if hasfield(typeof(cf), :lattice_angles)
        # Set α
        rotateAboutB!(view(primitive_vectors,3,:), 90u"°" - cf.lattice_angles[1])
        # Set β
        rotateAboutA!(view(primitive_vectors,3,:), 90u"°" - cf.lattice_angles[2])
        # Set γ
        rotateAboutC!(view(primitive_vectors,2,:), 90u"°" - cf.lattice_angles[3])
    end

    return primitive_vectors
end

FaceCenteredSupportedTypes = Union{Cubic, Orthorhombic}
function get_primitive_vectors(cf::FaceCenteredSupportedTypes, ct::FaceCentered)
    primitive_vectors = MMatrix{3,3}([0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0]) 
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)
    return primitive_vectors
end

BodyCenteredSupportedTypes = Union{Cubic, Orthorhombic, Tetragonal}
function get_primitive_vectors(cf::BodyCenteredSupportedTypes, ct::BodyCentered)
    primitive_vectors = MMatrix{3,3}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.5  0.5  0.5])
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)
    return primitive_vectors
end

BaseCenteredSupportedTypes = Union{Monoclinic, Orthorhombic}
function get_primitive_vectors(cf::BaseCenteredSupportedTypes, ct::BaseCentered)
    primitive_vectors = MMatrix{3,3}([1.0 1.0 0.0; 1.0 -1.0 0.0; 0.0  0.0  1.0])
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)
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

function get_primitive_vectors(cf::Rectangular, ct::Centered)
    primitive_vectors = MMatrix{2,2}([0.5 0.5; 0.5 -0.5])
    primitive_vectors = primitive_vectors.*transpose(cf.lattice_constants)
    return primitive_vectors
end


#####################################################
# Helper functions

rotateAboutA!(v, θ) = copyto!(v, MMatrix{3,3}([1.0 0.0 0.0; 0.0 cos(θ) -sin(θ); 0  sin(θ)  cos(θ)]) * v)
rotateAboutB!(v, θ) = copyto!(v, MMatrix{3,3}([cos(θ) 0.0 sin(θ); 0.0 1.0 0.0; -sin(θ) 0.0 cos(θ)]) * v)
rotateAboutC!(v, θ) = copyto!(v, MMatrix{3,3}([cos(θ) -sin(θ) 0.0; sin(θ) cos(θ) 0.0; 0.0  0.0  1.0]) * v)