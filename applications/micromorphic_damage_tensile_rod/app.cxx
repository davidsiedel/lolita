#include "lolita_lolita/lolita_core/lolita_core_n_4004.hxx"

int
main(int argc, char** argv)
{
    std::cout << std::fixed << std::setprecision(3);
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // declaring behavior
    auto lib_displacement_path = "/home/dsiedel/projetcs/lolita/lolita/applications/micromorphic_damage_tensile_rod/data/behavior/bhv_micromorphic_displacement/src/libBehaviour.so";
    auto lib_displacement_label = "MicromorphicDisplacement";
    // declaring behavior
    auto lib_damage_path = "/home/dsiedel/projetcs/lolita/lolita/applications/micromorphic_damage_tensile_rod/data/behavior/bhv_micromorphic_damage/src/libBehaviour.so";
    auto lib_damage_name = "MicromorphicDamage";
    //
    auto opts = mgis::behaviour::FiniteStrainBehaviourOptions{
        mgis::behaviour::FiniteStrainBehaviourOptions::PK1,
        mgis::behaviour::FiniteStrainBehaviourOptions::DPK1_DF
    };
    // auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // constants
    auto constexpr domain = lolita::Domain::cartesian(2);
    auto constexpr quadrature = lolita::Quadrature::gauss(2);
    auto constexpr cells = lolita::ElementType::cells(domain);
    auto constexpr faces = lolita::ElementType::faces(domain);
    // discretization
    auto constexpr hdg = lolita::HybridDiscontinuousGalerkin::hybridDiscontinuousGalerkin(1, 1);
    // generalized strains
    auto constexpr displacement_generalized_strain = lolita::GeneralizedStrain(lolita::Field::vector(), lolita::Mapping::smallStrain());
    auto constexpr damage_generalized_strain = lolita::GeneralizedStrain(lolita::Field::vector(), lolita::Mapping::smallStrain());
    // behaviors
    auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    auto constexpr damage_behavior = lolita::Behavior(damage_generalized_strain);
    // finite elements
    auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    auto constexpr damage_element =  lolita::FiniteElementMethod(damage_generalized_strain, damage_behavior, hdg, quadrature);
    // mesh
    auto file_path = "/home/dsiedel/projetcs/lolita/lolita/applications/micromorphic_damage_tensile_rod/mesh.msh";
    auto out_file = "/home/dsiedel/projetcs/lolita/lolita/applications/micromorphic_damage_tensile_rod/out.msh";
    // mesh build
    auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    // dofs
    auto face_displacement = elements->setDegreeOfFreedom<faces, lolita::Field::vector(), hdg.getFaceBasis()>("SQUARE", "Displacement");
    auto cell_displacement = elements->setDegreeOfFreedom<cells, lolita::Field::vector(), hdg.getCellBasis()>("SQUARE", "Displacement");
    auto face_damage = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("SQUARE", "Damage");
    auto cell_damage = elements->setDegreeOfFreedom<cells, lolita::Field::scalar(), hdg.getCellBasis()>("SQUARE", "Damage");
    //
    auto top_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
    auto left_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce");
    auto bottom_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
    // systems
    auto displacement_system = lolita::System::make();
    auto damage_system = lolita::System::make();
    displacement_system->setUnknown("Displacement", face_displacement->size());
    displacement_system->setBinding("TopForce", top_force->size());
    displacement_system->setBinding("LeftForce", left_force->size());
    displacement_system->setBinding("BottomForce", bottom_force->size());
    displacement_system->initialize();
    damage_system->setUnknown("Damage", face_damage->size());
    damage_system->initialize();
    // load
    auto load0 = elements->setConstraint<faces>("TOP", "Pull", [](lolita::Point const & p, lolita::Real const & t) { return t; }, 1, 0);
    auto load1 = elements->setConstraint<faces>("LEFT", "FixedL", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
    auto load2 = elements->setConstraint<faces>("BOTTOM", "FixedB", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 1, 0);
    // adding behavior
    auto micromorphic_displacement = elements->setBehavior<cells, quadrature>("SQUARE", lib_displacement_path, lib_displacement_label, hyp);
    auto micromorphic_damage = elements->setBehavior<cells, quadrature>("SQUARE", lib_damage_path, lib_damage_name, hyp);
    //making operators
    elements->setStrainOperators<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement");
    elements->setStrainOperators<cells, damage_element, hdg>("SQUARE", "MicromorphicDamage", "Damage");
    elements->setElementOperators<cells, displacement_element, hdg>("SQUARE", "DisplacementStabilization");
    elements->setElementOperators<cells, damage_element, hdg>("SQUARE", "DamageStabilization");
    // setting variable
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDisplacement", "YoungModulus", [](lolita::Point const & p) { return 1.0; });
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDisplacement", "PoissonRatio", [](lolita::Point const & p) { return 0.2; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDisplacement", "Temperature", [](lolita::Point const & p) { return 293.15; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDisplacement", "Damage", [](lolita::Point const & p) { return 0.0; });
    // setting variable
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDamage", "FractureEnergy", [](lolita::Point const & p) { return 1.0; });
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDamage", "CharacteristicLength", [](lolita::Point const & p) { return 0.1; });
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDamage", "PenalisationFactor", [](lolita::Point const & p) { return 300.; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDamage", "Temperature", [](lolita::Point const & p) { return 293.15; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDamage", "EnergyReleaseRate", [](lolita::Point const & p) { return 0.0; });
    // setting parameter
    elements->setParameter<faces>("TOP", "Lagrange", [](lolita::Point const & p) { return 1.0; });
    elements->setParameter<faces>("LEFT", "Lagrange", [](lolita::Point const & p) { return 1.0; });
    elements->setParameter<faces>("BOTTOM", "Lagrange", [](lolita::Point const & p) { return 1.0; });
    // stab
    elements->setParameter<cells>("SQUARE", "DisplacementStabilization", [](lolita::Point const & p) { return 1.0 / (1.0 + 0.2); });
    elements->setParameter<cells>("SQUARE", "DamageStabilization", [](lolita::Point const & p) { return 1.0 / 0.1; });
    lolita::GmshFileParser::setOutput<domain>(out_file, elements, "MicromorphicDisplacement");
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------

    auto load_max = 0.01;
    auto times = std::vector<lolita::Real>();
    for (auto i = 0; i < 11; i++)
    {
        times.push_back(i * load_max / 10);
    }
    
    auto step = 0;
    auto time = times[step];
    //
    auto displacement_newton_step = [&] ()
    {
        auto iteration = 0;
        while(iteration < 10)
        {
            displacement_system->initializeRhs();
            displacement_system->initializeLhs();
            displacement_system->initializeNormalization();
            elements->setStrainValues<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement");
            std::cout << "displacement iteration : " << iteration << std::endl;
            auto res_eval = elements->integrate<cells>("SQUARE", "MicromorphicDisplacement");
            if (res_eval.isFailure())
            {
                std::cout << "displacement integration failure" << std::endl;
                return false;
            }
            else
            {
                elements->assembleUnknownVector<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement", displacement_system);
                elements->assembleBindingVector<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system, time);
                elements->assembleBindingVector<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system, time);
                elements->assembleBindingVector<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system, time);
                std::cout << "displacement res eval : " << displacement_system->getResidualEvaluation() << std::endl;
                auto res = displacement_system->getResidualEvaluation();
                // std::cout << "res : " << lolita::Matrix<lolita::Real, 1, -1>(displacement_system->rhs_values_) << std::endl;
                if (res < 1.e-6)
                {
                    std::cout << "displacement convergence" << std::endl;
                    return true;
                }
                else
                {
                    std::cout << "displacement solve" << std::endl;
                    elements->assembleUnknownBlock<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement", displacement_system);
                    //
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system);
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system);
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system);
                    displacement_system->setCorrection();
                    elements->updateUnknown<cells, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
                    elements->updateUnknown<faces, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce", displacement_system);
            lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "MicromorphicDisplacement", 0);
            lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "MicromorphicDisplacement", 1);
            lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, step, time, "Displacement", "MicromorphicDisplacement", 0, 0);
            lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, step, time, "Displacement", "MicromorphicDisplacement", 1, 0);
                }
            }
            iteration ++;
        }
        std::cout << "displacement max iters" << std::endl;
        return false;
    };
    //
    auto damage_newton_step = [&] ()
    {
        auto iteration = 0;
        while(iteration < 10)
        {
            damage_system->initializeRhs();
            damage_system->initializeLhs();
            damage_system->initializeNormalization(1);
            elements->setStrainValues<cells, damage_element, hdg>("SQUARE", "MicromorphicDamage", "Damage");
            std::cout << "damage iteration : " << iteration << std::endl;
            auto res_eval = elements->integrate<cells>("SQUARE", "MicromorphicDamage");
            if (res_eval.isFailure())
            {
                std::cout << "damage integration failure" << std::endl;
                return false;
            }
            else
            {
                elements->assembleUnknownVector<cells, damage_element, hdg>("SQUARE", "MicromorphicDamage", "Damage", damage_system);
                // elements->assembleBindingVector<faces, damage_element, hdg>("TOP", "TopForce", "Damage", "Pull", displacement_system, time);
                // elements->assembleBindingVector<faces, damage_element, hdg>("LEFT", "LeftForce", "Damage", "FixedL", displacement_system, time);
                // elements->assembleBindingVector<faces, damage_element, hdg>("BOTTOM", "BottomForce", "Damage", "FixedB", displacement_system, time);
                std::cout << "damage res eval : " << damage_system->getResidualEvaluation() << std::endl;
                auto res = damage_system->getResidualEvaluation();
                // std::cout << "res : " << lolita::Matrix<lolita::Real, 1, -1>(displacement_system->rhs_values_) << std::endl;
                if (res < 1.e-6)
                {
                    std::cout << "damage convergence" << std::endl;
                    return true;
                }
                else
                {
                    std::cout << "damage solve" << std::endl;
                    elements->assembleUnknownBlock<cells, damage_element, hdg>("SQUARE", "MicromorphicDamage", "Damage", damage_system);
                    //
                    // elements->assembleBindingBlock<faces, damage_element, hdg>("TOP", "TopForce", "Damage", "Pull", displacement_system);
                    // elements->assembleBindingBlock<faces, damage_element, hdg>("LEFT", "LeftForce", "Damage", "FixedL", displacement_system);
                    // elements->assembleBindingBlock<faces, damage_element, hdg>("BOTTOM", "BottomForce", "Damage", "FixedB", displacement_system);
                    damage_system->setCorrection();
                    elements->updateUnknown<cells, damage_element, hdg>("SQUARE", "Damage", damage_system);
                    elements->updateUnknown<faces, damage_element, hdg>("SQUARE", "Damage", damage_system);
                    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce", displacement_system);
                    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce", displacement_system);
                    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce", displacement_system);
            lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "MicromorphicDamage", 0);
            lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "MicromorphicDamage", 1);
            lolita::GmshFileParser::addQuadratureDofOutput<2, domain, damage_element, hdg>(out_file, elements, step, time, "Damage", "MicromorphicDamage", 0, 0);
                }
            }
            iteration ++;
        }
        std::cout << "damage max iters" << std::endl;
        return false;
    };
    //
    auto displacement_convergence = [&] ()
    {
        displacement_system->initializeRhs();
        displacement_system->initializeNormalization();
        auto res_eval = elements->integrate<cells>("SQUARE", "MicromorphicDisplacement");
        if (res_eval.isFailure())
        {
            std::cout << "displacement integration failure" << std::endl;
            return -1;
        }
        else
        {
            elements->assembleUnknownVector<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement", displacement_system);
            elements->assembleBindingVector<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system, time);
            elements->assembleBindingVector<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system, time);
            elements->assembleBindingVector<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system, time);
            std::cout << "displacement res eval : " << displacement_system->getResidualEvaluation() << std::endl;
            auto res = displacement_system->getResidualEvaluation();
            // std::cout << "res : " << lolita::Matrix<lolita::Real, 1, -1>(displacement_system->rhs_values_) << std::endl;
            if (res < 1.e-6)
            {
                std::cout << "displacement convergence" << std::endl;
                return 1;
            }
            else
            {
                return 0;
            }
        }
    };
    //
    auto coupled_newton_step = [&] ()
    {
        auto glob_iter = 0;
        while(glob_iter < 1000)
        {
            if (displacement_newton_step())
            {
                auto set_energy = [&] (
                    auto const & finite_element
                )
                {
                    auto cnt = 0;
                    for (auto & ip : finite_element->quadrature_.at("MicromorphicDamage").ips_)
                    {
                        auto const & dip = finite_element->quadrature_.at("MicromorphicDisplacement").ips_[cnt];
                        auto const * rhs = mgis::behaviour::getInternalStateVariable(dip.behavior_data_->s1, "EnergyReleaseRate");
                        auto * lhs = mgis::behaviour::getExternalStateVariable(ip.behavior_data_->s1, "EnergyReleaseRate");
                        * lhs = * rhs;
                        // mgis::behaviour::setExternalStateVariable(ip.behavior_data_->s1, "EnergyReleaseRate", * rhs);
                        cnt ++;
                    }
                };
                elements->caller<cells>("SQUARE", set_energy);
            }
            else
            {
                return false;
            }
            if (damage_newton_step())
            {
                auto set_energy = [&] (
                    auto const & finite_element
                )
                {
                    auto cnt = 0;
                    for (auto & ip : finite_element->quadrature_.at("MicromorphicDisplacement").ips_)
                    {
                        auto const & dip = finite_element->quadrature_.at("MicromorphicDamage").ips_[cnt];
                        auto const * rhs = mgis::behaviour::getInternalStateVariable(dip.behavior_data_->s1, "Damage");
                        auto * lhs = mgis::behaviour::getExternalStateVariable(ip.behavior_data_->s1, "Damage");
                        * lhs = * rhs;
                        // mgis::behaviour::setExternalStateVariable(ip.behavior_data_->s1, "Damage", * rhs);
                        cnt ++;
                    }
                };
                elements->caller<cells>("SQUARE", set_energy);
            }
            else
            {
                return false;
            }
            auto conv_test = displacement_convergence();
            if (conv_test == 1)
            {
                return true;
            }
            else if (conv_test == -1)
            {
                return false;
            }
            glob_iter ++;
        }
        return false;
    };
    //
    while (step < times.size())
    {
        if (coupled_newton_step())
        {
            std::cout << "time step convergence" << std::endl;
            elements->reserveBehaviorData<cells>("SQUARE", "MicromorphicDisplacement");
            elements->reserveBehaviorData<cells>("SQUARE", "MicromorphicDamage");
            //
            elements->reserveUnknownCoefficients<cells, lolita::Field::vector(), hdg.getCellBasis()>("SQUARE", "Displacement");
            elements->reserveUnknownCoefficients<faces, lolita::Field::vector(), hdg.getFaceBasis()>("SQUARE", "Displacement");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
            //
            elements->reserveUnknownCoefficients<cells, lolita::Field::scalar(), hdg.getCellBasis()>("SQUARE", "Damage");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("SQUARE", "Damage");
            // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "Elasticity", 0);
            // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_file, elements, step, time, "Elasticity", 1);
            // lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, step, time, "Displacement", "Elasticity", 0, 0);
            // lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_file, elements, step, time, "Displacement", "Elasticity", 1, 0);
            step ++;
            time = times[step];
            /* code */
        }
        else
        {
            std::cout << "time step split" << std::endl;
            time = times[step - 1] + (1.0 / 2.0) * (time - times[step - 1]);
        }
    }
}