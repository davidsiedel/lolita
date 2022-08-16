#include "lolita_lolita/lolita_core/lolita_core_n_5000.hxx"

int main(int argc, char** argv)
{
    // constants
    // auto constexpr domain = lolita::Domain::cartesian(2);
    // auto constexpr cell_basis = lolita::Basis::monomial(1);
    // auto constexpr face_basis = lolita::Basis::monomial(1);
    // auto constexpr quadrature = lolita::Quadrature::gauss(2);
    // auto constexpr cells = lolita::ElementType::cells(domain);
    // auto constexpr faces = lolita::ElementType::faces(domain);
    // // fields
    // auto constexpr displacement_field = lolita::Field::vector();
    // auto constexpr damage_field = lolita::Field::scalar();
    // // generalized strains
    // auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(displacement_field, lolita::Mapping::smallStrain());
    // auto constexpr damage_generalized_strain = lolita::GeneralizedStrain(damage_field, lolita::Mapping::gradient(), lolita::Mapping::identity());
    // // behaviors
    // auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    // auto constexpr damage_behavior = lolita::Behavior(damage_generalized_strain);
    // // discretization
    // auto constexpr hdg = lolita::HybridDiscontinuousGalerkin(cell_basis, face_basis, lolita::HybridDiscontinuousGalerkin::Stabilization::Hdg);
    // // finite elements
    // auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    // auto constexpr damage_element =  lolita::FiniteElementMethod(damage_generalized_strain, damage_behavior, hdg, quadrature);
    // // mesh    
    // auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/data/meshes/unit_square_3_cpp.msh";
    // // dofs
    // auto face_displacement = std::make_shared<lolita::DegreeOfFreedom>("FaceDisplacement");
    // auto cell_displacement = std::make_shared<lolita::DegreeOfFreedom>("CellDisplacement");
    // auto face_damage = std::make_shared<lolita::DegreeOfFreedom>("FaceDamage");
    // auto cell_damage = std::make_shared<lolita::DegreeOfFreedom>("CellDamage");
    // // lag
    // auto top_force = std::make_shared<lolita::DegreeOfFreedom>("TopForce");
    // // load
    // auto load_f = std::make_shared<lolita::Load>([](lolita::Point const &p, lolita::Real const &t) { return 1.0; }, 0, 0);
    // // mesh build
    // auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    // elements->addElement<cells>("Displacement", "SQUARE");
    // elements->addElement<faces>("Displacement", "SQUARE");
    // elements->addElement<cells>("Damage", "SQUARE");
    // elements->addElement<faces>("Damage", "SQUARE");
    // // dofs
    // elements->addDegreeOfFreedom<cells, displacement_field, face_basis>("Displacement", "SQUARE", cell_displacement);
    // elements->addDegreeOfFreedom<faces, displacement_field, face_basis>("Displacement", "SQUARE", face_displacement);
    // elements->addDegreeOfFreedom<cells, damage_field, face_basis>("Damage", "SQUARE", cell_damage);
    // elements->addDegreeOfFreedom<faces, damage_field, face_basis>("Damage", "SQUARE", face_damage);
    // // load
    // elements->addLoad<faces>("Displacement", "TOP", load_f);
    // elements->addDegreeOfFreedom<faces, displacement_field, face_basis>("Displacement", "TOP", top_force);
    // // bhv
    // auto lib_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/behaviour/src/libBehaviour.so";
    // auto lib_name = "Voce";
    // auto opts = mgis::behaviour::FiniteStrainBehaviourOptions{
    //     mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
    //     mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF
    // };
    // auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;

    // auto bhvv = std::make_shared<mgis::behaviour::Behaviour>(mgis::behaviour::load(opts, lib_path, lib_name, hyp));
    // elements->addBehavior<cells, quadrature>("SQUARE", bhvv);
    // elements->addBehavior<cells>("Displacement", "SQUARE", lib_name);

}