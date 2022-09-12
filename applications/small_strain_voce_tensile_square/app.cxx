#include "core/700_mesh.hxx"

int
main(int argc, char** argv)
{
    std::cout << "Have " << argc << " arguments:" << std::endl;
    for (int i = 0; i < argc; ++i) {
        std::cout << argv[i] << std::endl;
    }
    std::cout << std::fixed << std::setprecision(3);
    std::cout << "hardware_concurrency : " << std::thread::hardware_concurrency() << std::endl;
    // auto time_steps = lolita::utility::linspace(0., 1., 10);
    // for (auto val : time_steps)
    // {
    //     std::cout << val << std::endl;
    // }
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // declaring behavior
    auto lib_path = "/home/dsiedel/projetcs/lolita/applications/linear_elastic_tensile_square/data/behavior/bhv_elasticity/src/libBehaviour.so";
    // auto lib_path = "data/behavior/bhv_elasticity/src/libBehaviour.so";
    auto lib_name = "Elasticity";
    auto opts = mgis::behaviour::FiniteStrainBehaviourOptions{
        mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
        mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF
    };
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    // constants
    auto constexpr domain = lolita::Domain::cartesian(2);
    auto constexpr quadrature = lolita::Quadrature::gauss(2);
    auto constexpr cells = lolita::ElementType::cells(domain);
    auto constexpr faces = lolita::ElementType::faces(domain);
    // fields
    auto constexpr displacement_field = lolita::Field::vector();
    // generalized strains
    auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(lolita::Field::vector(), lolita::Mapping::smallStrain());
    // behaviors
    auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    // discretization
    auto constexpr hdg = lolita::HybridDiscontinuousGalerkin::hybridDiscontinuousGalerkin(1, 1);
    // finite elements
    auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    // mesh
    auto file_path = "/home/dsiedel/projetcs/lolita/applications/small_strain_voce_tensile_square/mesh2.msh";
    auto out_file = "/home/dsiedel/projetcs/lolita/applications/small_strain_voce_tensile_square/out.msh";
    // auto file_path = "mesh.msh";
    // auto out_file = "out.msh";
    // mesh build
    auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    // -> DEBUG
    // elements->getElements<2, 1>()[1]->template getSymmetricGradientRhsDEBUG<lolita::Field::vector(), hdg>(0, 0);
    // <- DEBUG
    // dofs
    auto face_displacement = elements->setDegreeOfFreedom<faces, lolita::Field::vector(), hdg.getFaceBasis()>("SQUARE", "Displacement");
    auto cell_displacement = elements->setDegreeOfFreedom<cells, lolita::Field::vector(), hdg.getCellBasis()>("SQUARE", "Displacement");
    //
    auto top_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
    auto left_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce");
    auto bottom_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
    // systems
    auto displacement_system = lolita::System::make();
    displacement_system->setUnknown("Displacement", face_displacement->size());
    displacement_system->setBinding("TopForce", top_force->size());
    displacement_system->setBinding("LeftForce", left_force->size());
    displacement_system->setBinding("BottomForce", bottom_force->size());
    displacement_system->initialize();
    std::cout << "U size : " << displacement_system->getUnknownsSize() << std::endl;
    std::cout << "B size : " << displacement_system->getBindingsSize() << std::endl;
    // load
    auto load0 = elements->setConstraint<faces>("TOP", "Pull", [](lolita::Point const & p, lolita::Real const & t) { return t; }, 1, 0);
    auto load1 = elements->setConstraint<faces>("LEFT", "FixedL", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
    auto load2 = elements->setConstraint<faces>("BOTTOM", "FixedB", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 1, 0);
    // adding behavior
    auto micromorphic_damage = elements->setBehavior<cells, quadrature>("SQUARE", lib_path, lib_name, hyp);
    //making operators
    elements->setStrainOperators<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement");
    elements->setElementOperators<cells, displacement_element, hdg>("SQUARE", "DisplacementStabilization");
    // setting variable
    // elements->setMaterialProperty<cells>("SQUARE", "Elasticity", "YoungModulus", [](lolita::Point const & p) { return 1.0; });
    // elements->setMaterialProperty<cells>("SQUARE", "Elasticity", "PoissonRatio", [](lolita::Point const & p) { return 0.0; });
    elements->setExternalVariable<cells>("SQUARE", "Elasticity", "Temperature", [](lolita::Point const & p) { return 293.15; });
    // elements->setExternalVariable<cells>("SQUARE", "Elasticity", "Damage", [](lolita::Point const & p) { return 0.0; });
    // setting parameter
    elements->setParameter<faces>("TOP", "Lagrange", [](lolita::Point const & p) { return 206.9e9; });
    elements->setParameter<faces>("LEFT", "Lagrange", [](lolita::Point const & p) { return 206.9e9; });
    elements->setParameter<faces>("BOTTOM", "Lagrange", [](lolita::Point const & p) { return 206.9e9; });
    // stab
    elements->setParameter<cells>("SQUARE", "DisplacementStabilization", [](lolita::Point const & p) { return 206.9e9; });
    lolita::GmshFileParser::setOutput<domain>(out_file, elements, "Elasticity");
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------

    // -> DEBUG
    // auto elem_count = 0;
    // for (auto const & finite_element : elements->getElements<2, 0>())
    // {
    //     std::cout << "****** element :" << std::endl;
    //     std::cout << lolita::mat2str(finite_element->getCurrentCentroid()) << std::endl;
    //     std::cout << "-- stabilization :" << std::endl;
    //     std::cout << lolita::mat2str(finite_element->operators_.at("DisplacementStabilization")) << std::endl;
    //     auto ip_count = 0;
    //     for (auto const & ip : finite_element->quadrature_.at("Elasticity").ips_)
    //     {
    //         auto grad = finite_element->template getMapping<lolita::Field::vector(), lolita::Mapping::smallStrain(), hdg>(ip.reference_coordinates_);
    //         std::cout << "-- grad " << ip_count << std::endl;
    //         std::cout << lolita::mat2str(ip.ops_.at("Displacement")) << std::endl;
    //         ip_count ++;
    //     }
    //     elem_count ++;
    // }
    // for (auto const & finite_element : elements->getElements<2, 1>())
    // {
    //     std::cout << "****** element :" << std::endl;
    //     std::cout << lolita::mat2str(finite_element->getCurrentCentroid()) << std::endl;
    //     std::cout << "-- stabilization :" << std::endl;
    //     std::cout << lolita::mat2str(finite_element->operators_.at("DisplacementStabilization")) << std::endl;
    //     auto ip_count = 0;
    //     for (auto const & ip : finite_element->quadrature_.at("Elasticity").ips_)
    //     {
    //         auto grad = finite_element->template getMapping<lolita::Field::vector(), lolita::Mapping::smallStrain(), hdg>(ip.reference_coordinates_);
    //         std::cout << "-- grad " << ip_count << std::endl;
    //         std::cout << lolita::mat2str(ip.ops_.at("Displacement")) << std::endl;
    //         ip_count ++;
    //     }
    //     elem_count ++;
    // }
    // <- DEBUG

    auto num_steps = 5;
    auto times = std::vector<lolita::Real>();
    for (auto i = 0; i < 5 + 1; i++)
    {
        times.push_back(i * 1.0e-2 / 5);
    }
    
    auto step = 0;
    auto time = times[step];
    auto newton_step = [&] ()
    {
        auto iteration = 0;
        while(iteration < 10)
        {
            displacement_system->initialize();
            elements->setStrainValues<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement");
            std::cout << "iteration : " << iteration << std::endl;
            auto res_eval = elements->integrate<cells>("SQUARE", "Elasticity");
            if (res_eval.isFailure())
            {
                std::cout << "integration failure" << std::endl;
                return false;
            }
            else
            {
                elements->assembleUnknownVector<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement", displacement_system);
                elements->assembleBindingVector<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system, time);
                elements->assembleBindingVector<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system, time);
                elements->assembleBindingVector<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system, time);
                std::cout << "res eval : " << displacement_system->getResidualEvaluation() << std::endl;
                auto res = displacement_system->getResidualEvaluation();
                // std::cout << "rhs : " << lolita::mat2str(displacement_system->rhs_values_) << std::endl;
                if (res < 1.e-6)
                {
                    std::cout << "convergence" << std::endl;
                    return true;
                }
                else
                {
                    elements->assembleUnknownBlock<cells, displacement_element, hdg>("SQUARE", "Elasticity", "Displacement", displacement_system);
                    //
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system);
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system);
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system);
                    // std::cout << "solve" << std::endl;
                    displacement_system->setCorrection();
                    // std::cout << lolita::mat2str(displacement_system->correction_values_) << std::endl;
                    elements->updateUnknown<cells, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
                    elements->updateUnknown<faces, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce", displacement_system);
            // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "Elasticity", 0);
            // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "Elasticity", 1);
            // lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, step, time, "Displacement", "Elasticity", 0, 0);
            // lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, step, time, "Displacement", "Elasticity", 1, 0);
                }
            }
            iteration ++;
        }
        std::cout << "max iters" << std::endl;
        return false;
    };
    while (step < times.size())
    {
        // -> DEBUG
        // break;
        // <- DEBUG
        std::cout << "step : " << step << " time : " << time << std::endl;
        if (newton_step())
        {
            std::cout << "time step convergence" << std::endl;
            elements->reserveBehaviorData<cells>("SQUARE", "Elasticity");
            elements->reserveUnknownCoefficients<cells, lolita::Field::vector(), hdg.getCellBasis()>("SQUARE", "Displacement");
            elements->reserveUnknownCoefficients<faces, lolita::Field::vector(), hdg.getFaceBasis()>("SQUARE", "Displacement");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
            lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "Elasticity", 0);
            lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "Elasticity", 1);
            lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, step, time, "Displacement", "Elasticity", 0, 0);
            lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, step, time, "Displacement", "Elasticity", 1, 0);
            step ++;
            time = times[step];
            /* code */
        }
        else
        {
            std::cout << "time step split" << std::endl;
            time = times[step - 1] + (1.0 / 2.0) * (time - times[step - 1]);
            // -> DEBUG
            // step ++;
            // time = times[step];
            // <- DEBUG
        }
    }
}