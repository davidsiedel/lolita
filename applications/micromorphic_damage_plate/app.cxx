#include "lolita_lolita/lolita_core/lolita_core_n_4004.hxx"

// #include "petscsys.h" /* framework routines */
// #include "petscvec.h" /* vectors */
// #include "petscmat.h" /* matrices */
// #include "petscksp.h"
// #include <mpi.h>

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
    auto lib_displacement_path = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_plate/data/behavior/bhv_micromorphic_displacement/src/libBehaviour.so";
    auto lib_displacement_label = "MicromorphicDisplacement";
    // declaring behavior
    auto lib_damage_path = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_plate/data/behavior/bhv_micromorphic_damage/src/libBehaviour.so";
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
    auto file_path = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_plate/mesh.msh";
    auto out_displacement_file = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_plate/out_u.msh";
    auto out_damage_file = "/home/dsiedel/projetcs/lolita/applications/micromorphic_damage_plate/out_d.msh";
    // mesh build
    auto elements = lolita::MeshFileParser(file_path).template makeFiniteElementSet<domain>();
    // dofs
    auto face_displacement = elements->setDegreeOfFreedom<faces, lolita::Field::vector(), hdg.getFaceBasis()>("SQUARE", "Displacement");
    auto cell_displacement = elements->setDegreeOfFreedom<cells, lolita::Field::vector(), hdg.getCellBasis()>("SQUARE", "Displacement");
    auto face_damage = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("SQUARE", "Damage");
    auto cell_damage = elements->setDegreeOfFreedom<cells, lolita::Field::scalar(), hdg.getCellBasis()>("SQUARE", "Damage");
    //
    auto top_force = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
    auto circle_force_x = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceX");
    auto circle_force_y = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceY");
    auto circle_force_d = elements->setDegreeOfFreedom<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceD");
    // systems
    auto displacement_system = lolita::System::make();
    displacement_system->setUnknown("Displacement", face_displacement->size());
    displacement_system->setBinding("TopForce", top_force->size());
    displacement_system->setBinding("CircleForceX", circle_force_x->size());
    displacement_system->setBinding("CircleForceY", circle_force_y->size());
    auto damage_system = lolita::System::make();
    damage_system->setUnknown("Damage", face_damage->size());
    damage_system->setBinding("CircleForceD", circle_force_d->size());
    // load
    auto load0 = elements->setConstraint<faces>("TOP", "Pull", [](lolita::Point const & p, lolita::Real const & t) { return t; }, 1, 0);
    auto load1 = elements->setConstraint<faces>("CIRCLE", "FixedX", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
    auto load2 = elements->setConstraint<faces>("CIRCLE", "FixedY", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 1, 0);
    auto load3 = elements->setConstraint<faces>("CIRCLE", "FixedD", [](lolita::Point const & p, lolita::Real const & t) { return 0.0; }, 0, 0);
    // adding behavior
    auto micromorphic_displacement = elements->setBehavior<cells, quadrature>("SQUARE", lib_displacement_path, lib_displacement_label, hyp);
    auto micromorphic_damage = elements->setBehavior<cells, quadrature>("SQUARE", lib_damage_path, lib_damage_name, hyp);
    //making operators
    elements->setStrainOperators<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement");
    elements->setStrainOperators<cells, damage_element, hdg>("SQUARE", "MicromorphicDamage", "Damage");
    elements->setElementOperators<cells, displacement_element, hdg>("SQUARE", "DisplacementStabilization");
    elements->setElementOperators<cells, damage_element, hdg>("SQUARE", "DamageStabilization");
    // setting variable
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDisplacement", "YoungModulus", [](lolita::Point const & p) { return 200.0; });
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDisplacement", "PoissonRatio", [](lolita::Point const & p) { return 0.2; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDisplacement", "Temperature", [](lolita::Point const & p) { return 293.15; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDisplacement", "Damage", [](lolita::Point const & p) { return 0.0; });
    // setting variable
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDamage", "FractureEnergy", [](lolita::Point const & p) { return 1.0; });
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDamage", "CharacteristicLength", [](lolita::Point const & p) { return 0.02; });
    elements->setMaterialProperty<cells>("SQUARE", "MicromorphicDamage", "PenalisationFactor", [](lolita::Point const & p) { return 300.; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDamage", "Temperature", [](lolita::Point const & p) { return 293.15; });
    elements->setExternalVariable<cells>("SQUARE", "MicromorphicDamage", "EnergyReleaseRate", [](lolita::Point const & p) { return 0.0; });
    // setting parameter
    elements->setParameter<faces>("TOP", "TopForceLagrange", [](lolita::Point const & p) { return 200.0; });
    elements->setParameter<faces>("CIRCLE", "CircleForceXLagrange", [](lolita::Point const & p) { return 200.0; });
    elements->setParameter<faces>("CIRCLE", "CircleForceYLagrange", [](lolita::Point const & p) { return 200.0; });
    elements->setParameter<faces>("CIRCLE", "CircleForceDLagrange", [](lolita::Point const & p) { return 1.0; });
    // stab
    elements->setParameter<cells>("SQUARE", "DisplacementStabilization", [](lolita::Point const & p) { return 200.0 / (1.0 + 0.2); });
    elements->setParameter<cells>("SQUARE", "DamageStabilization", [](lolita::Point const & p) { return 1.0 / 0.05; });
    //
    lolita::GmshFileParser::setOutput<domain>(out_displacement_file, elements, "MicromorphicDisplacement");
    lolita::GmshFileParser::setOutput<domain>(out_damage_file, elements, "MicromorphicDamage");
    //
    elements->setBandWidth<faces, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
    elements->setBandWidth<faces, damage_element, hdg>("SQUARE", "Displacement", damage_system);
    std::cout << "setBandWidth : " << displacement_system->band_width_ << std::endl;
    // -> TEST
    // std::cout << "elements->getElements<2, 1>()[0]->quadrature_.at(MicromorphicDamage).ips_[0].ops_.at(Damage)" << std::endl;
    // std::cout << elements->getElements<2, 1>()[0]->quadrature_.at("MicromorphicDamage").ips_[0].ops_.at("Damage") << std::endl;
    // std::cout << "elements->getElements<2, 1>()[0]->quadrature_.at(MicromorphicDisplacement).ips_[0].ops_.at(Displacement)" << std::endl;
    // std::cout << elements->getElements<2, 1>()[0]->quadrature_.at("MicromorphicDisplacement").ips_[0].ops_.at("Displacement") << std::endl;
    // <- TEST
    
    // auto num_steps = 50;
    // std::cout << "times :" << std::endl;
    // auto times = std::vector<lolita::Real>();
    // for (auto i = 0; i < num_steps + 1; i++)
    // {
    //     auto val = i * 2.0e-1 / num_steps;
    //     std::cout << i << " : " << val << " / ";
    //     times.push_back(val);
    // }
    // std::cout << std::endl;

    auto times = std::vector<lolita::Real>{
        0.0, 0.014000000000000002, 0.028000000000000004, 0.04200000000000001, 0.05600000000000001, 0.07, 0.0722, 0.07440000000000001, 0.0766, 0.07880000000000001, 0.081, 0.08320000000000001, 0.0854, 0.08760000000000001, 0.0898, 0.092, 0.0942, 0.0964, 0.09860000000000001, 0.1008, 0.10300000000000001, 0.1052, 0.1074, 0.1096, 0.11180000000000001, 0.114, 0.1162, 0.1184, 0.1206, 0.12279999999999999, 0.125
    };

    auto tick = std::chrono::high_resolution_clock::now();
    auto tock = std::chrono::high_resolution_clock::now();
    
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
            tick = std::chrono::high_resolution_clock::now();
            elements->setStrainValues<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement");
            tock = std::chrono::high_resolution_clock::now();
            auto strain_time = std::chrono::duration<double>(tock - tick);
            tick = std::chrono::high_resolution_clock::now();
            auto res_eval = elements->integrate<cells>("SQUARE", "MicromorphicDisplacement");
            tock = std::chrono::high_resolution_clock::now();
            auto integration_time = std::chrono::duration<double>(tock - tick);
            // std::cout << "---------------------------------------------> INTEGRATION" << std::endl;
            // std::cout << "--- strain    u : " << std::setprecision(10) << std::scientific << strain_time.count() << std::endl;
            // std::cout << "--- integrate u : " << std::setprecision(10) << std::scientific << integration_time.count() << std::endl;
            // std::cout << "---------------------------------------------< INTEGRATION" << std::endl;
            if (res_eval.isFailure())
            {
                std::cout << "!!! displacement integration failure" << std::endl;
                return false;
            }
            else
            {
                tick = std::chrono::high_resolution_clock::now();
                elements->assembleUnknownVector<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement", displacement_system);
                tock = std::chrono::high_resolution_clock::now();
                auto u_v_assembly_time = std::chrono::duration<double>(tock - tick);
                tick = std::chrono::high_resolution_clock::now();
                elements->assembleBindingVector<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system, time);
                elements->assembleBindingVector<faces, displacement_element, hdg>("CIRCLE", "CircleForceX", "Displacement", "FixedX", displacement_system, time);
                elements->assembleBindingVector<faces, displacement_element, hdg>("CIRCLE", "CircleForceY", "Displacement", "FixedY", displacement_system, time);
                tock = std::chrono::high_resolution_clock::now();
                auto b_v_assembly_time = std::chrono::duration<double>(tock - tick);
                // std::cout << "---------------------------------------------> RHS" << std::endl;
                // std::cout << "--- RHS u : " << std::setprecision(10) << std::scientific << u_v_assembly_time.count() << std::endl;
                // std::cout << "--- RHS b : " << std::setprecision(10) << std::scientific << b_v_assembly_time.count() << std::endl;
                // std::cout << "---------------------------------------------< RHS" << std::endl;
                auto res = displacement_system->getResidualEvaluation();
                if (res < 1.e-6)
                {
                    std::cout << "displacement iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | convergence" << std::endl;
                    return true;
                }
                else
                {
                    std::cout << "displacement iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | solve" << std::endl;
                    tick = std::chrono::high_resolution_clock::now();
                    elements->assembleUnknownBlock<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement", displacement_system);
                    tock = std::chrono::high_resolution_clock::now();
                    auto u_assembly_time = std::chrono::duration<double>(tock - tick);
                    //
                    tick = std::chrono::high_resolution_clock::now();
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system);
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("CIRCLE", "CircleForceX", "Displacement", "FixedX", displacement_system);
                    elements->assembleBindingBlock<faces, displacement_element, hdg>("CIRCLE", "CircleForceY", "Displacement", "FixedY", displacement_system);
                    tock = std::chrono::high_resolution_clock::now();
                    auto b_assembly_time = std::chrono::duration<double>(tock - tick);
                    // std::cout << "solve" << std::endl;
                    tick = std::chrono::high_resolution_clock::now();
                    displacement_system->setCorrection();
                    tock = std::chrono::high_resolution_clock::now();
                    auto u_solver_time = std::chrono::duration<double>(tock - tick);
                    // std::cout << lolita::mat2str(displacement_system->correction_values_) << std::endl;
                    tick = std::chrono::high_resolution_clock::now();
                    elements->updateUnknown<cells, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
                    elements->updateUnknown<faces, displacement_element, hdg>("SQUARE", "Displacement", displacement_system);
                    // * face_displacement += displacement_system->getUnknownCorrection("Displacement");
                    tock = std::chrono::high_resolution_clock::now();
                    auto u_update_time = std::chrono::duration<double>(tock - tick);
                    tick = std::chrono::high_resolution_clock::now();
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceX", displacement_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceY", displacement_system);
                    tock = std::chrono::high_resolution_clock::now();
                    auto b_update_time = std::chrono::duration<double>(tock - tick);
                    // std::cout << "---------------------------------------------> SOLVER" << std::endl;
                    // std::cout << "--- LHS u : " << std::setprecision(10) << std::scientific << u_assembly_time.count() << std::endl;
                    // std::cout << "--- LHS b : " << std::setprecision(10) << std::scientific << b_assembly_time.count() << std::endl;
                    // std::cout << "--- solver   u : " << std::setprecision(10) << std::scientific << u_solver_time.count() << std::endl;
                    // std::cout << "--- update   u : " << std::setprecision(10) << std::scientific << u_update_time.count() << std::endl;
                    // std::cout << "--- update   b : " << std::setprecision(10) << std::scientific << b_update_time.count() << std::endl;
                    // std::cout << "---------------------------------------------< SOLVER" << std::endl;
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
            elements->setStrainValues<cells, damage_element, hdg>("SQUARE", "MicromorphicDamage", "Damage");
            // std::cout << "damage iteration : " << iteration << std::endl;
            auto res_eval = elements->integrate<cells>("SQUARE", "MicromorphicDamage");
            if (res_eval.isFailure())
            {
                std::cout << "!!! damage integration failure" << std::endl;
                return false;
            }
            else
            {
                elements->assembleUnknownVector<cells, damage_element, hdg>("SQUARE", "MicromorphicDamage", "Damage", damage_system);
                elements->assembleBindingVector<faces, damage_element, hdg>("CIRCLE", "CircleForceD", "Damage", "FixedD", damage_system, time);
                auto res = damage_system->getResidualEvaluation();
                if (res < 1.e-6)
                {
                    auto val = elements->getDissipatedEnergy<cells>("SQUARE", "MicromorphicDamage");
                    std::cout << "---------------------------------------------> DISSIPATED ENERGY" << std::endl;
                    std::cout << "--- : " << std::setprecision(10) << std::scientific << val << std::endl;
                    std::cout << "---------------------------------------------< DISSIPATED ENERGY" << std::endl;
                    std::cout << "damage iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | convergence" << std::endl;
                    // std::cout << "damage convergence" << std::endl;
                    return true;
                }
                else
                {
                    std::cout << "damage iteration " << iteration << " | res " << std::setprecision(10) << std::scientific << res << " | solve" << std::endl;
                    elements->assembleUnknownBlock<cells, damage_element, hdg>("SQUARE", "MicromorphicDamage", "Damage", damage_system);
                    elements->assembleBindingBlock<faces, damage_element, hdg>("CIRCLE", "CircleForceD", "Damage", "FixedD", damage_system);
                    // std::cout << "solve" << std::endl;
                    damage_system->setCorrection();
                    // std::cout << lolita::mat2str(damage_system->correction_values_) << std::endl;
                    elements->updateUnknown<cells, damage_element, hdg>("SQUARE", "Damage", damage_system);
                    elements->updateUnknown<faces, damage_element, hdg>("SQUARE", "Damage", damage_system);
                    elements->updateBinding<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceD", damage_system);
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
        auto res_eval = elements->integrate<cells>("SQUARE", "MicromorphicDisplacement");
        if (res_eval.isFailure())
        {
            std::cout << "!!! displacement integration failure" << std::endl;
            return -1;
        }
        else
        {
            elements->assembleUnknownVector<cells, displacement_element, hdg>("SQUARE", "MicromorphicDisplacement", "Displacement", displacement_system);
            elements->assembleBindingVector<faces, displacement_element, hdg>("TOP", "TopForce", "Displacement", "Pull", displacement_system, time);
            elements->assembleBindingVector<faces, displacement_element, hdg>("CIRCLE", "CircleForceX", "Displacement", "FixedX", displacement_system, time);
            elements->assembleBindingVector<faces, displacement_element, hdg>("CIRCLE", "CircleForceY", "Displacement", "FixedY", displacement_system, time);
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
        // -> DEBUG
        // break;
        // <- DEBUG
        std::cout << "**** step : " << step << " time : " << time << std::endl;
        if (coupled_newton_step())
        {
            std::cout << "-- time step convergence" << std::endl;
            elements->reserveBehaviorData<cells>("SQUARE", "MicromorphicDisplacement");
            elements->reserveUnknownCoefficients<cells, lolita::Field::vector(), hdg.getCellBasis()>("SQUARE", "Displacement");
            elements->reserveUnknownCoefficients<faces, lolita::Field::vector(), hdg.getFaceBasis()>("SQUARE", "Displacement");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceX");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceY");
            //
            elements->reserveBehaviorData<cells>("SQUARE", "MicromorphicDamage");
            elements->reserveUnknownCoefficients<cells, lolita::Field::scalar(), hdg.getCellBasis()>("SQUARE", "Damage");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("SQUARE", "Damage");
            elements->reserveUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceD");
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
                lolita::GmshFileParser::addQuadratureDofOutput<2, domain, displacement_element, hdg>(out_damage_file, elements, step, time, "Damage", "MicromorphicDamage", 0, 0);
            #endif
            step ++;
            time = times[step];
            /* code */
        }
        else
        {
            std::cout << "-- time step split" << std::endl;
            elements->recoverBehaviorData<cells>("SQUARE", "MicromorphicDisplacement");
            elements->recoverUnknownCoefficients<cells, lolita::Field::vector(), hdg.getCellBasis()>("SQUARE", "Displacement");
            elements->recoverUnknownCoefficients<faces, lolita::Field::vector(), hdg.getFaceBasis()>("SQUARE", "Displacement");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("TOP", "TopForce");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceX");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceY");
            //
            elements->recoverBehaviorData<cells>("SQUARE", "MicromorphicDamage");
            elements->recoverUnknownCoefficients<cells, lolita::Field::scalar(), hdg.getCellBasis()>("SQUARE", "Damage");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("SQUARE", "Damage");
            elements->recoverUnknownCoefficients<faces, lolita::Field::scalar(), hdg.getFaceBasis()>("CIRCLE", "CircleForceD");
            time = times[step - 1] + (1.0 / 2.0) * (time - times[step - 1]);
            // -> DEBUG
            // step ++;
            // time = times[step];
            // <- DEBUG
        }
    }
}