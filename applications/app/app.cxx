#include "lolita_lolita/lolita_core/lolita_core_n_4004.hxx"

int main(int argc, char** argv)
{

    // tick = std::chrono::high_resolution_clock::now();

    // // declaring behavior
    // auto lib_path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/bhv_elasticity/src/libBehaviour.so";
    // auto lib_name = "Elasticity";
    // auto opts = mgis::behaviour::FiniteStrainBehaviourOptions{
    //     mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
    //     mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF
    // };
    // auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    // //
    // //
    // //
    // // constants
    // auto constexpr domain = lolita::Domain::cartesian(2);
    // auto constexpr quadrature = lolita::Quadrature::gauss(2);
    // auto constexpr cells = lolita::ElementType::cells(domain);
    // auto constexpr faces = lolita::ElementType::faces(domain);
    // // fields
    // auto constexpr displacement_field = lolita::Field::vector();
    // // generalized strains
    // auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(lolita::Field::vector(), lolita::Mapping::smallStrain());
    // // behaviors
    // auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    // // discretization
    // auto constexpr hdg = lolita::HybridDiscontinuousGalerkin::hybridDiscontinuousGalerkin(1, 1);
    // // finite elements
    // auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    // // mesh
    // auto file_path = "/home/dsiedel/projetcs/lolita/lolita/tests/input/mesh.msh";
    // //
    // auto out_file = "/home/dsiedel/projetcs/lolita/lolita/tests/output/out.msh";
    // // mesh build
    // auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    // // dofs
    // auto face_displacement = elements->setDegreeOfFreedom<faces, lolita::Field::vector(), hdg.getFaceBasis()>("SQUARE", "Displacement");
    // auto cell_displacement = elements->setDegreeOfFreedom<cells, lolita::Field::vector(), hdg.getCellBasis()>("SQUARE", "Displacement");
    // //
    // auto top_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
    // auto left_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce");
    // auto bottom_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
    // // systems
    // auto displacement_system = lolita::System::make();
    // displacement_system->setUnknown("Displacement", face_displacement->size());
    // displacement_system->setBinding("TopForce", top_force->size());
    // displacement_system->setBinding("LeftForce", left_force->size());
    // displacement_system->setBinding("BottomForce", bottom_force->size());
    // displacement_system->initialize();
    // std::cout << "displacement_system u size : " << displacement_system->getUnknownsSize() << std::endl;
    // std::cout << "displacement_system b size : " << displacement_system->getBindingsSize() << std::endl;
    // // load
    // auto load0 = elements->setConstraint<faces>("TOP", "Pull", [](lolita::Point const & p, lolita::Real const & t) { return 1.0; }, 1, 0);
    // auto load1 = elements->setConstraint<faces>("LEFT", "FixedL", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
    // auto load2 = elements->setConstraint<faces>("BOTTOM", "FixedB", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 1, 0);
    // // adding behavior
    // auto micromorphic_damage = elements->setBehavior<cells, quadrature>("SQUARE", lib_path, lib_name, hyp);
    // //making operators
    // elements->setStrainOperators<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement");
    // elements->setElementOperators<cells, displacement_element, hdg>("SQUARE", "Stabilization");
    // // setting variable
    // elements->setExternalVariable<cells>("SQUARE", "Elasticity", "Temperature", [](lolita::Point const & p) { return 293.15; });
    // // setting parameter
    // elements->setParameter<faces>("TOP", "Lagrange", [](lolita::Point const & p) { return 1.0; });
    // elements->setParameter<faces>("LEFT", "Lagrange", [](lolita::Point const & p) { return 1.0; });
    // elements->setParameter<faces>("BOTTOM", "Lagrange", [](lolita::Point const & p) { return 1.0; });
    // // stab
    // elements->setParameter<cells>("SQUARE", "Stabilization", [](lolita::Point const & p) { return 1.0; });
    // //
    // elements->setStrainValues<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement");
    // //
    // elements->integrate<cells>("SQUARE", "Elasticity");
    // //
    // elements->assembleUnknownBlock<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement", displacement_system);
    // //
    // elements->assembleBindingBlock<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system);
    // elements->assembleBindingBlock<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system);
    // elements->assembleBindingBlock<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system);
    // //
    // displacement_system->setCorrection();
    // //
    // elements->updateUnknown<cells, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
    // elements->updateUnknown<faces, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce", displacement_system);
    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce", displacement_system);
    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce", displacement_system);
    // //
    // lolita::GmshFileParser::setOutput<domain>(out_file, elements, "Elasticity");
    // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, 0, 0.0, "Elasticity", 0);
    // lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, 0, 0.0, "Displacement", "Elasticity", 0, 0);
    // lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, 0, 0.0, "Displacement", "Elasticity", 1, 0);

}