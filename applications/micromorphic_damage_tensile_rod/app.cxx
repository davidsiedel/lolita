#include "lolita_lolita/lolita_core/lolita_core_n_4004.hxx"

// class Polygon {
//     public:
//     virtual void show() {
//         std::cout<<"It is a polygon"<<std::endl;
//     }
// };

// class Hexagon : public Polygon {
//     public:
//     void show() {
//         std::cout<<"Hexagon is a 6 sided polygon"<<std::endl;
//     }
// };

// class Pentagon : public Polygon {
//     public:
//     void show() {
//         std::cout<<"Pentagon is a 5 sided polygon"<<std::endl;
//     }
// };

int
main(int argc, char** argv)
{
    // std::cout << std::fixed << std::setprecision(3);
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // declaring behavior
    auto lib_displacement_path = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/data/behavior/bhv_micromorphic_displacement/src/libBehaviour.so";
    auto lib_displacement_label = "MicromorphicDisplacement";
    // declaring behavior
    auto lib_damage_path = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/data/behavior/bhv_micromorphic_damage/src/libBehaviour.so";
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
    auto constexpr damage_generalized_strain = lolita::GeneralizedStrain(lolita::Field::scalar(), lolita::Mapping::gradient(), lolita::Mapping::identity());
    // behaviors
    auto constexpr displacement_behavior = lolita::Behavior(displacement_generalized_strain);
    auto constexpr damage_behavior = lolita::Behavior(damage_generalized_strain);
    // finite elements
    auto constexpr displacement_element =  lolita::FiniteElementMethod(displacement_generalized_strain, displacement_behavior, hdg, quadrature);
    auto constexpr damage_element =  lolita::FiniteElementMethod(damage_generalized_strain, damage_behavior, hdg, quadrature);
    // mesh
    auto file_path = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/mesh.msh";
    auto out_displacement_file = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/out_u.msh";
    auto out_damage_file = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/out_d.msh";
    auto force_out_stream = std::basic_ofstream<lolita::Character>();
    auto damage_dissipated_energy_out_stream = std::basic_ofstream<lolita::Character>();
    auto damage_stored_energy_out_stream = std::basic_ofstream<lolita::Character>();
    auto displacement_dissipated_energy_out_stream = std::basic_ofstream<lolita::Character>();
    auto displacement_stored_energy_out_stream = std::basic_ofstream<lolita::Character>();
    force_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/out_force.txt");
    damage_dissipated_energy_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/out_damage_dissipated_energy.txt");
    damage_stored_energy_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/out_damage_stored_energy.txt");
    displacement_dissipated_energy_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/out_displacement_dissipated_energy.txt");
    displacement_stored_energy_out_stream.open("/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_tensile_rod/out_displacement_stored_energy.txt");
    // mesh build
    auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    // dofs
    auto face_displacement = elements->setDegreeOfFreedom<faces, lolita::Field::vector(), hdg.getFaceBasis()>("ROD", "Displacement");
    auto cell_displacement = elements->setDegreeOfFreedom<cells, lolita::Field::vector(), hdg.getCellBasis()>("ROD", "Displacement");
    auto face_damage = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("ROD", "Damage");
    auto cell_damage = elements->setDegreeOfFreedom<cells, lolita::Field::scalar(), hdg.getCellBasis()>("ROD", "Damage");
    //
    auto top_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
    auto left_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce");
    // <- RIGHT
    auto right_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("RIGHT", "RightForce");
    // <- RIGHT
    auto bottom_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
    // systems
    auto displacement_system = lolita::System::make();
    auto damage_system = lolita::System::make();
    displacement_system->setUnknown("Displacement", face_displacement->size());
    displacement_system->setBinding("TopForce", top_force->size());
    displacement_system->setBinding("LeftForce", left_force->size());
    //
    displacement_system->setBinding("RightForce", right_force->size());
    //
    displacement_system->setBinding("BottomForce", bottom_force->size());
    displacement_system->initialize();
    damage_system->setUnknown("Damage", face_damage->size());
    damage_system->initialize();
    // load
    auto load0 = elements->setConstraint<faces>("TOP", "Pull", [](lolita::Point const & p, lolita::Real const & t) { return t; }, 1, 0);
    auto load1 = elements->setConstraint<faces>("LEFT", "FixedL", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
    // -> RIGHT
    auto loadr = elements->setConstraint<faces>("RIGHT", "FixedR", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
    // <- RIGHT
    auto load2 = elements->setConstraint<faces>("BOTTOM", "FixedB", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 1, 0);
    // adding behavior
    auto micromorphic_displacement = elements->setBehavior<cells, quadrature>("ROD", lib_displacement_path, lib_displacement_label, hyp);
    auto micromorphic_damage = elements->setBehavior<cells, quadrature>("ROD", lib_damage_path, lib_damage_name, hyp);
    //making operators
    elements->setStrainOperators<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement");
    elements->setStrainOperators<cells, damage_element, hdg>("ROD", "MicromorphicDamage", "Damage");
    elements->setElementOperators<cells, displacement_element, hdg>("ROD", "DisplacementStabilization");
    elements->setElementOperators<cells, damage_element, hdg>("ROD", "DamageStabilization");
    // setting variable
    elements->setMaterialProperty<cells>("SANE", "MicromorphicDisplacement", "YoungModulus", [](lolita::Point const & p) { return 200.0; });
    elements->setMaterialProperty<cells>("DEFECT", "MicromorphicDisplacement", "YoungModulus", [](lolita::Point const & p) { return 199.5; });
    elements->setMaterialProperty<cells>("ROD", "MicromorphicDisplacement", "PoissonRatio", [](lolita::Point const & p) { return 0.0; });
    elements->setExternalVariable<cells>("ROD", "MicromorphicDisplacement", "Temperature", [](lolita::Point const & p) { return 293.15; });
    elements->setExternalVariable<cells>("ROD", "MicromorphicDisplacement", "Damage", [](lolita::Point const & p) { return 0.0; });
    // setting variable
    elements->setMaterialProperty<cells>("ROD", "MicromorphicDamage", "FractureEnergy", [](lolita::Point const & p) { return 1.0; });
    elements->setMaterialProperty<cells>("ROD", "MicromorphicDamage", "CharacteristicLength", [](lolita::Point const & p) { return 0.05; });
    elements->setMaterialProperty<cells>("ROD", "MicromorphicDamage", "PenalisationFactor", [](lolita::Point const & p) { return 300.; });
    elements->setExternalVariable<cells>("ROD", "MicromorphicDamage", "Temperature", [](lolita::Point const & p) { return 293.15; });
    elements->setExternalVariable<cells>("ROD", "MicromorphicDamage", "EnergyReleaseRate", [](lolita::Point const & p) { return 0.0; });
    // setting parameter
    elements->setParameter<faces>("TOP", "TopForceLagrange", [](lolita::Point const & p) { return 200.0; });
    elements->setParameter<faces>("LEFT", "LeftForceLagrange", [](lolita::Point const & p) { return 200.0; });
    // <- RIGHT
    elements->setParameter<faces>("RIGHT", "RightForceLagrange", [](lolita::Point const & p) { return 200.0; });
    // <- RIGHT
    elements->setParameter<faces>("BOTTOM", "BottomForceLagrange", [](lolita::Point const & p) { return 200.0; });
    // stab
    elements->setParameter<cells>("ROD", "DisplacementStabilization", [](lolita::Point const & p) { return 200.0 / (1.0 + 0.0); });
    elements->setParameter<cells>("ROD", "DamageStabilization", [](lolita::Point const & p) { return 1.0 / 0.05; });
    //
    lolita::GmshFileParser::setOutput<domain>(out_displacement_file, elements, "MicromorphicDisplacement");
    lolita::GmshFileParser::setOutput<domain>(out_damage_file, elements, "MicromorphicDamage");
    // -> TEST
    std::cout << "elements->getElements<2, 1>()[0]->quadrature_.at(MicromorphicDamage).ips_[0].ops_.at(Damage)" << std::endl;
    std::cout << elements->getElements<2, 1>()[0]->quadrature_.at("MicromorphicDamage").ips_[0].ops_.at("Damage") << std::endl;
    std::cout << "elements->getElements<2, 1>()[0]->quadrature_.at(MicromorphicDisplacement).ips_[0].ops_.at(Displacement)" << std::endl;
    std::cout << elements->getElements<2, 1>()[0]->quadrature_.at("MicromorphicDisplacement").ips_[0].ops_.at("Displacement") << std::endl;
    // <- TEST
    
    auto num_steps = 50;
    std::cout << "times :" << std::endl;
    auto times = std::vector<lolita::Real>();
    for (auto i = 0; i < num_steps + 1; i++)
    {
        auto val = i * 3.0e-1 / num_steps;
        std::cout << i << " : " << val << " / ";
        times.push_back(val);
    }
    std::cout << std::endl;
    
    auto step = 0;
    auto time = times[step];
    // -> DISPLACEMENT
    auto displacement_newton_step = [&] ()
    {
        auto iteration = 0;
        while(iteration < 10)
        {
            // displacement_system->initialize();
            displacement_system->initializeRhs();
            displacement_system->initializeLhs();
            displacement_system->initializeNormalization();
            elements->setStrainValues<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement");
            auto res_eval = elements->integrate<cells>("ROD", "MicromorphicDisplacement");
            if (res_eval.isFailure())
            {
                std::cout << "!!! displacement integration failure" << std::endl;
                return false;
            }
            else
            {
                elements->assembleUnknownVector<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement", displacement_system);
                elements->assembleBindingVector<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system, time);
                elements->assembleBindingVector<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system, time);
                // <- RIGHT
                elements->assembleBindingVector<faces, displacement_element, hdg>("RIGHT", "RightForce", "Displacement", "FixedR", displacement_system, time);
                // -> RIGHT
                elements->assembleBindingVector<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system, time);
                // std::cout << "displacement res eval : " << std::setprecision(17) << displacement_system->getResidualEvaluation() << std::endl;
                // std::cout << "displacement res norm : " << std::setprecision(17) << displacement_system->getNormalization() << std::endl;
                auto res = displacement_system->getResidualEvaluation();
                // std::cout << "rhs\n" << lolita::mat2str(displacement_system->rhs_values_) << std::endl;
                if (res < 1.e-6)
                {
                    std::cout << "displacement iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | convergence" << std::endl;
                    return true;
                }
                else
                {
                    std::cout << "displacement iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | solve" << std::endl;
                    elements->assembleUnknownBlock<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement", displacement_system);
                    //
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system);
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system);
                    // <- RIGHT
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("RIGHT", "RightForce", "Displacement", "FixedR", displacement_system);
                    // -> RIGHT
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system);
                    // std::cout << "solve" << std::endl;
                    displacement_system->setCorrection();
                    // std::cout << lolita::mat2str(displacement_system->correction_values_) << std::endl;
                    elements->updateUnknown<cells, displacement_element, hdg>("ROD", "Displacement", displacement_system);
                    elements->updateUnknown<faces, displacement_element, hdg>("ROD", "Displacement", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce", displacement_system);
                    // <- RIGHT
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("RIGHT", "RightForce", displacement_system);
                    // -> RIGHT
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce", displacement_system);
                    #ifdef DEBUG
                        // std::cout << "writing debug output" << std::endl;
                        // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_displacement_file, elements, total_iteration, total_iteration, "MicromorphicDisplacement", 0);
                        // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_displacement_file, elements, total_iteration, total_iteration, "MicromorphicDisplacement", 1);
                        // lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_displacement_file, elements, total_iteration, total_iteration, "MicromorphicDisplacement", 0);
                        // lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_displacement_file, elements, total_iteration, total_iteration, "MicromorphicDisplacement", 1);
                        // lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_displacement_file, elements, total_iteration, total_iteration, "Displacement", "MicromorphicDisplacement", 0, 0);
                        // lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_displacement_file, elements, total_iteration, total_iteration, "Displacement", "MicromorphicDisplacement", 1, 0);
                    #endif
                }
            }
            iteration ++;
            // total_iteration ++;
        }
        std::cout << "!!! displacement max iters" << std::endl;
        return false;
    };
    // <- DISPLACEMENT
    // -> DAMAGE
    auto damage_newton_step = [&] ()
    {
        auto iteration = 0;
        while(iteration < 10)
        {
            // displacement_system->initialize();
            damage_system->initializeRhs();
            damage_system->initializeLhs();
            damage_system->initializeNormalization(1);
            elements->setStrainValues<cells, damage_element, hdg>("ROD", "MicromorphicDamage", "Damage");
            // std::cout << "damage iteration : " << iteration << std::endl;
            auto res_eval = elements->integrate<cells>("ROD", "MicromorphicDamage");
            if (res_eval.isFailure())
            {
                std::cout << "!!! damage integration failure" << std::endl;
                return false;
            }
            else
            {
                elements->assembleUnknownVector<cells, damage_element, hdg>("ROD", "MicromorphicDamage", "Damage", damage_system);
                // elements->assembleBindingVector<faces, damage_element, hdg>("TOP", "TopForce", "Damage", "Pull", damage_system, time);
                // elements->assembleBindingVector<faces, damage_element, hdg>("LEFT", "LeftForce", "Damage", "FixedL", damage_system, time);
                // // <- RIGHT
                // elements->assembleBindingVector<faces, damage_element, hdg>("RIGHT", "RightForce", "Damage", "FixedR", damage_system, time);
                // // -> RIGHT
                // elements->assembleBindingVector<faces, damage_element, hdg>("BOTTOM", "BottomForce", "Damage", "FixedB", damage_system, time);
                // std::cout << "damage res eval : " << std::setprecision(17) << damage_system->getResidualEvaluation() << std::endl;
                // std::cout << "damage res norm : " << std::setprecision(17) << damage_system->getNormalization() << std::endl;
                auto res = damage_system->getResidualEvaluation();
                // std::cout << "displacement iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | convergence" << std::endl;
                // std::cout << "rhs\n" << lolita::mat2str(damage_system->rhs_values_) << std::endl;
                if (res < 1.e-6)
                {
                    std::cout << "damage iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | convergence" << std::endl;
                    // std::cout << "damage convergence" << std::endl;
                    return true;
                }
                else
                {
                    std::cout << "damage iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | solve" << std::endl;
                    elements->assembleUnknownBlock<cells, damage_element, hdg>("ROD", "MicromorphicDamage", "Damage", damage_system);
                    //
                    // elements->assembleBindingBlock<faces, damage_element, hdg>("TOP", "TopForce", "Damage", "Pull", damage_system);
                    // elements->assembleBindingBlock<faces, damage_element, hdg>("LEFT", "LeftForce", "Damage", "FixedL", damage_system);
                    // // <- RIGHT
                    // elements->assembleBindingBlock<faces, damage_element, hdg>("RIGHT", "RightForce", "Damage", "FixedR", damage_system);
                    // // -> RIGHT
                    // elements->assembleBindingBlock<faces, damage_element, hdg>("BOTTOM", "BottomForce", "Damage", "FixedB", damage_system);
                    // std::cout << "solve" << std::endl;
                    damage_system->setCorrection();
                    // std::cout << lolita::mat2str(damage_system->correction_values_) << std::endl;
                    elements->updateUnknown<cells, damage_element, hdg>("ROD", "Damage", damage_system);
                    elements->updateUnknown<faces, damage_element, hdg>("ROD", "Damage", damage_system);
                    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce", damage_system);
                    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce", damage_system);
                    // // <- RIGHT
                    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("RIGHT", "RightForce", damage_system);
                    // // -> RIGHT
                    // elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce", damage_system);
                    #ifdef DEBUG
                        // std::cout << "writing debug output" << std::endl;
                        // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_damage_file, elements, total_iteration, total_iteration, "MicromorphicDamage", 0);
                        // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_damage_file, elements, total_iteration, total_iteration, "MicromorphicDamage", 1);
                        // lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_damage_file, elements, total_iteration, total_iteration, "MicromorphicDamage", 2);
                        // lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_damage_file, elements, total_iteration, total_iteration, "MicromorphicDamage", 0);
                        // lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_damage_file, elements, total_iteration, total_iteration, "MicromorphicDamage", 1);
                        // lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_damage_file, elements, total_iteration, total_iteration, "MicromorphicDamage", 2);
                        // lolita::GmshFileParser::addQuadratureDofOutput<2, domain, damage_element, hdg>(out_damage_file, elements, total_iteration, total_iteration, "Damage", "MicromorphicDamage", 0, 0);
                    #endif
                }
            }
            iteration ++;
            // total_iteration ++;
        }
        std::cout << "!!! damage max iters" << std::endl;
        return false;
    };
    // <- DAMAGE
    auto displacement_convergence = [&] ()
    {
        displacement_system->initializeRhs();
        displacement_system->initializeNormalization();
        auto res_eval = elements->integrate<cells>("ROD", "MicromorphicDisplacement");
        if (res_eval.isFailure())
        {
            std::cout << "!!! displacement integration failure" << std::endl;
            return -1;
        }
        else
        {
            elements->assembleUnknownVector<cells, displacement_element, hdg>("ROD", "MicromorphicDisplacement", "Displacement", displacement_system);
            elements->assembleBindingVector<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system, time);
            elements->assembleBindingVector<faces, displacement_element, hdg>("LEFT", "LeftForce", "Displacement", "FixedL", displacement_system, time);
            elements->assembleBindingVector<faces, displacement_element, hdg>("RIGHT", "RightForce", "Displacement", "FixedR", displacement_system, time);
            elements->assembleBindingVector<faces, displacement_element, hdg>("BOTTOM", "BottomForce", "Displacement", "FixedB", displacement_system, time);
            // std::cout << "displacement res eval : " << displacement_system->getResidualEvaluation() << std::endl;
            auto res = displacement_system->getResidualEvaluation();
            // std::cout << "res : " << lolita::Matrix<lolita::Real, 1, -1>(displacement_system->rhs_values_) << std::endl;
            if (res < 1.e-6)
            {
                // std::cout << "step convergence" << std::endl;
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
                elements->caller<cells>("ROD", set_energy);
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
                elements->caller<cells>("ROD", set_energy);
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
        // -> DEBUG
        // break;
        // <- DEBUG
        std::cout << "**** step : " << step << " time : " << time << std::endl;
        if (coupled_newton_step())
        {
            std::cout << "-- time step convergence" << std::endl;
            auto top_force_value = elements->getBindingIntegral<faces, displacement_element, hdg>("TOP", "TopForce", 0, 0);
            auto damage_stored_energy_value = elements->getStoredEnergy<cells>("ROD", "MicromorphicDamage");
            auto damage_dissipated_energy_value = elements->getDissipatedEnergy<cells>("ROD", "MicromorphicDamage");
            auto displacement_stored_energy_value = elements->getStoredEnergy<cells>("ROD", "MicromorphicDisplacement");
            auto displacement_dissipated_energy_value = elements->getDissipatedEnergy<cells>("ROD", "MicromorphicDisplacement");
force_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << top_force_value << "\n";
damage_stored_energy_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << damage_stored_energy_value << "\n";
damage_dissipated_energy_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << damage_dissipated_energy_value << "\n";
displacement_stored_energy_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << displacement_stored_energy_value << "\n";
displacement_dissipated_energy_out_stream << step << ", "<< time << ", " << std::setprecision(10) << std::scientific << displacement_dissipated_energy_value << "\n";
            elements->reserveBehaviorData<cells>("ROD", "MicromorphicDisplacement");
            elements->reserveUnknownCoefficients<cells, lolita::Field::vector(), hdg.getCellBasis()>("ROD", "Displacement");
            elements->reserveUnknownCoefficients<faces, lolita::Field::vector(), hdg.getFaceBasis()>("ROD", "Displacement");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce");
            // <- RIGHT
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("RIGHT", "RightForce");
            // -> RIGHT
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
            //
            elements->reserveBehaviorData<cells>("ROD", "MicromorphicDamage");
            elements->reserveUnknownCoefficients<cells, lolita::Field::scalar(), hdg.getCellBasis()>("ROD", "Damage");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("ROD", "Damage");
            #ifndef DEBUG
                // std::cout << "writing regular output" << std::endl;
                lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_displacement_file, elements, step, time, "MicromorphicDisplacement", 0);
                lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_displacement_file, elements, step, time, "MicromorphicDisplacement", 1);
                lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_displacement_file, elements, step, time, "MicromorphicDisplacement", 0);
                lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_displacement_file, elements, step, time, "MicromorphicDisplacement", 1);
                lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_displacement_file, elements, step, time, "Displacement", "MicromorphicDisplacement", 0, 0);
                lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_displacement_file, elements, step, time, "Displacement", "MicromorphicDisplacement", 1, 0);
                //
                lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 0);
                lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 1);
                lolita::GmshFileParser::addQuadratureStrainOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 2);
                lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 0);
                lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 1);
                lolita::GmshFileParser::addQuadratureStressOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 2);
                lolita::GmshFileParser::addQuadratureInternalVariableOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 0);
                // lolita::GmshFileParser::addQuadratureInternalVariableOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 1);
                // lolita::GmshFileParser::addQuadratureInternalVariableOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 2);
                // lolita::GmshFileParser::addQuadratureInternalVariableOutput<2, domain>(out_damage_file, elements, step, time, "MicromorphicDamage", 3);
                lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_damage_file, elements, step, time, "Damage", "MicromorphicDamage", 0, 0);
            #endif
            step ++;
            time = times[step];
            /* code */
        }
        else
        {
            std::cout << "-- time step split" << std::endl;
            elements->recoverBehaviorData<cells>("ROD", "MicromorphicDisplacement");
            elements->recoverUnknownCoefficients<cells, lolita::Field::vector(), hdg.getCellBasis()>("ROD", "Displacement");
            elements->recoverUnknownCoefficients<faces, lolita::Field::vector(), hdg.getFaceBasis()>("ROD", "Displacement");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("LEFT", "LeftForce");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("RIGHT", "RightForce");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("BOTTOM", "BottomForce");
            //
            elements->recoverBehaviorData<cells>("ROD", "MicromorphicDamage");
            elements->recoverUnknownCoefficients<cells, lolita::Field::scalar(), hdg.getCellBasis()>("ROD", "Damage");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("ROD", "Damage");
            time = times[step - 1] + (1.0 / 2.0) * (time - times[step - 1]);
            // -> DEBUG
            // step ++;
            // time = times[step];
            // <- DEBUG
        }
    }
}