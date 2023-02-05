#Add support for vectors in N_unit_cells
 
#Only supports 3d crystals, cant make something like graphene


export
    Crystal

 
#########
# TYEPS #
#########

abstract type Lattice end
abstract type CrystalFamily end
abstract type CenteringType end

#Centering types for multiple dispatch
struct Primitive <: CenteringType end
struct FaceCentered <: CenteringType end
struct BodyCentered <: CenteringType end
struct BaseCentered <: CenteringType end

struct BravaisLattice{D} <: Lattice
    crystal_family::CrystalFamily
    centering_type::CenteringType
    primitive_vectors::MMatrix{D,D}
end

function BravaisLattice(cf::CrystalFamily, ct::CenteringType, dim::Integer)
    p_vec = get_primitive_vectors(cf,ct)
    return BravaisLattice{dim}(cf, ct, p_vec)
end


struct Atom{D}
    sym::Symbol
    position::SVector{D}
end

struct BasisAtom{D}
    basis_vector::SVector{D}
    atom::Atom
end


struct Crystal{D}
    lattice::BravaisLattice{D}
    basis::SVector{BasisAtom{D}}
end



#####################################################

struct Cubic{LC} <: CrystalFamily
    lattice_constants::SVector{3,LC}
end

function Cubic(a)
    return Cubic{typeof(lattice_constant)}(SVector(a,a,a))
end

#####################################################

struct Orthorhombic{LC} <: CrystalFamily
    lattice_constants::SVector{3,LC}
end

function Orthorhombic(a,b,c)
    return Orthorhombic{typeof(a)}(SVector(a,b,c))
end

#####################################################

struct Monoclinic{LC,LA} <: CrystalFamily
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

function Monoclinic(a, b, c, β)
    return Monoclinic{typeof(a),typeof(β)}(SVector(a,b,c), SVector(90u"°", β, 90u"°"))
end

#####################################################

struct Triclinic{LC,LA} <: CrystalFamily
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

function Triclinic(a, b, c, α, β, γ)
    return Triclinic{typeof(a),typeof(β)}(SVector(a,b,c), SVector(α, β, γ))
end
#####################################################

struct Tetragonal{LC} <: CrystalFamily
    lattice_constants::SVector{3,LC}
end

function Tetragonal(a, c)
    return Tetragonal{typeof(a)}(SVector(a,a,c))
end

#####################################################

struct Rhombohedral{LC,LA} <: CrystalFamily
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

function Rhombohedral(a, α)
    return Rhombohedral{typeof(a),typeof(α)}(SVector(a,a,a),SVector(α, α, α))
end

#####################################################

struct Hexagonal{LC,LA} <: CrystalFamily
    lattice_constants::SVector{3,LC}
    lattice_angles::SVector{3,LA}
end

function Hexagonal(a,c,γ)
    return HexagonalBravaisLattice{typeof(a),typeof(γ)}(SVector(a,a,c),SVector(90u"°", 90u"°", γ))
end

#####################################################

function get_primitive_vectors(cf::CrystalFamily, ct::Primitive)
    primitive_vectors = MMatrix{3,3}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0])
    primitive_vectors .*= transpose(cf.lattice_constants)

    if hasfield(cf, lattice_angles)
        #rotate a-axis
        rotateAboutB!(view(primitive_vectors,1,:), 90u"°" - cf.lattice_angles[2])
        rotateAboutC!(view(primitive_vectors,1,:), 90u"°" - cf.lattice_angles[3])
        #rotate b-axis
        rotateAboutA!(view(primitive_vectors,2,:) ,90u"°" - cf.lattice_angles[1])
        rotateAboutC!(view(primitive_vectors,2,:), 90u"°" - cf.lattice_angles[3])
        #rotate c-axis
        rotateAboutA!(view(primitive_vectors,3,:), 90u"°" - cf.lattice_angles[1])
        rotateAboutB!(view(primitive_vectors,3,:), 90u"°" - cf.lattice_angles[2])
    end

    return primitive_vectors
end

FaceCenteredSupportedTypes = Union{Cubic, Orthorhombic}
function get_primitive_vectors(cf::FaceCenteredSupportedTypes, ct::FaceCentered)
    primitive_vectors = MMatrix{3,3}([0.0 0.5 0.5; 0.5 0.0 0.5; 0.5 0.5 0.0])
    primitive_vectors .*= transpose(cf.lattice_constants)
    return primitive_vectors
end

BodyCenteredSupportedTypes = Union{Cubic, Orthorhombic, Tetragonal}
function get_primitive_vectors(cf::BodyCenteredSupportedTypes, ct::BodyCentered)
    primitive_vectors = MMatrix{3,3}([1.0 0.0 0.0; 0.0 1.0 0.0; 0.5  0.5  0.5])
    primitive_vectors .*= transpose(cf.lattice_constants)
    return primitive_vectors
end

BaseCenteredSupportedTypes = Union{Monoclinic, Orthorhombic}
function get_primitive_vectors(cf::BaseCenteredSupportedTypes, ct::BaseCentered)
    primitive_vectors = MMatrix{3,3}([1.0 1.0 0.0; 1.0 -1.0 0.0; 0.0  0.0  1.0])
    primitive_vectors .*= transpose(cf.lattice_constants)
    return primitive_vectors
end


#####################################################
# Helper functions

rotateAboutA!(v, θ) = copyto!(v, MMatrix{3,3}([1.0 0.0 0.0; 0.0 cos(θ) -sin(θ); 0  sin(θ)  cos(θ)]) * v)
rotateAboutB!(v, θ) = copyto!(v, MMatrix{3,3}([cos(θ) 0.0 sin(θ); 0.0 1.0 0.0; -sin(θ) 0.0 cos(θ)]) * v)
rotateAboutC!(v, θ) = copyto!(v, MMatrix{3,3}([cos(θ) -sin(θ) 0.0; sin(θ) cos(θ) 0.0; 0.0  0.0  1.0]) * v)
