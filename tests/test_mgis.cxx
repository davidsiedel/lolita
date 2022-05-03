//
// Created by dsiedel on 13/04/2022.
//

#include "gtest/gtest.h"

#include <span>

#include "lolita/lolita_containers.hxx"
#include "lolita/lolita.hxx"
#include "lolita/lolita_matrix.hxx"

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

TEST(test_mgis, test_mgis_1)
{

    using namespace lolita;

    auto path = "";
    path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/behaviour/src/libBehaviour.so";
    path = "/home/dsiedel/projetcs/lolita/lolita/tests/data/bhv_micromorphic/src/libBehaviour.so";
    auto name = "";
    name = "Voce";
    name = "MicromorphicDamageII";
    auto hyp = mgis::behaviour::Hypothesis::PLANESTRAIN;
    auto opt = mgis::behaviour::FiniteStrainBehaviourOptions{
        mgis::behaviour::FiniteStrainBehaviourOptions::StressMeasure::PK1,
        mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator::DPK1_DF
    };
//    opt.stress_measure = mgis::behaviour::FiniteStrainBehaviourOptions::PK1;
//    opt.tangent_operator = mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator::DPK1_DF;

//    auto bhv = mgis::behaviour::load(opt, path, name, hyp);
    auto bhv = mgis::behaviour::load(path, name, hyp);

    mgis::behaviour::InitialStateView s0;
    mgis::behaviour::StateView s1;

    mgis::behaviour::BehaviourDataView bdv;

//    mgis_bv_make_behaviour_data_view(&bdv,)
    auto bd = mgis::behaviour::BehaviourData(bhv);
    auto bdp = std::make_shared<mgis::behaviour::BehaviourData>(mgis::behaviour::BehaviourData(bhv));
    auto value = Real(245.0);
    mgis::behaviour::setMaterialProperty(bdp->s0, "FractureEnergy", value);
    mgis::behaviour::setMaterialProperty(bdp->s1, "FractureEnergy", value);
    mgis::behaviour::setMaterialProperty(bdp->s0, "CharacteristicLength", value);
    mgis::behaviour::setMaterialProperty(bdp->s1, "CharacteristicLength", value);
    mgis::behaviour::setMaterialProperty(bdp->s0, "PenalisationFactor", value);
    mgis::behaviour::setMaterialProperty(bdp->s1, "PenalisationFactor", value);
    mgis::behaviour::setExternalStateVariable(bdp->s0, "Temperature", value);
    mgis::behaviour::setExternalStateVariable(bdp->s1, "Temperature", value);
    mgis::behaviour::setExternalStateVariable(bdp->s0, "EnergyReleaseRate", value);
    mgis::behaviour::setExternalStateVariable(bdp->s1, "EnergyReleaseRate", value);

    auto start = std::chrono::high_resolution_clock::now();
    auto bdv2 = mgis::behaviour::make_view(* bdp);
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = duration_cast<std::chrono::milliseconds>(stop - start);
    bdp->K[0] = 4;
    Intg res0 = mgis::behaviour::integrate(bdv2, bhv);
    print("time (s) :", duration.count() * 1.e-3);

//    matrix::MatMap<Vector<Real, 3>>(bdp->s0.gradients.data()) = Vector<Real, 3>::Zero();
//    matrix::MatMap<Vector<Real, 3>>(bdp->s1.gradients.data()) = Vector<Real, 3>::Zero();

    print("K");

    for (int i = 0; i < bdp->K.size(); ++i) {
        print(bdp->K[i]);
    }

    print("grad");

    for (int i = 0; i < bdp->s0.gradients.size(); ++i) {
        print(bdp->s0.gradients[i]);
    }

    print("---");


    mgis::behaviour::MaterialDataManager material_data_manager(bhv, 2);
//    bdv.

    auto vals = std::vector<Real>{245.0, 245.0};
    auto vals2 = std::array<Real, 2>{245.0, 245.0};
    auto vals3 = std::tuple<Real, Real>{245.0, 245.0};
    auto const constexpr vals4 = std::array<Real, 2>{245.0, 245.0};
    auto valsspac = std::span(vals2);

    auto cstpr = [] <auto K> () constexpr
    {
        auto const constexpr vals4 = std::array<Real, 1>{K};
        auto valsspan = std::span(vals4);
        return Indx(valsspan[0]);
    };

    std::span<Real> myspan;

    auto const constexpr vlmp = cstpr.operator ()<2>();

    std::array<Real, vlmp> t;

    mgis::behaviour::setExternalStateVariable(
            material_data_manager.s0,
            "Temperature",
            valsspac
    );
    mgis::behaviour::setExternalStateVariable(
            material_data_manager.s1,
            "Temperature",
            245.0
    );

    mgis::behaviour::setExternalStateVariable(
            material_data_manager.s0,
            "EnergyReleaseRate",
            245.0
    );
    mgis::behaviour::setExternalStateVariable(
            material_data_manager.s1,
            "EnergyReleaseRate",
            245.0
    );

    mgis::behaviour::setMaterialProperty(
            material_data_manager.s0,
            "FractureEnergy",
            245.0
    );

    mgis::behaviour::setMaterialProperty(
            material_data_manager.s1,
            "FractureEnergy",
            245.0
    );

    mgis::behaviour::setMaterialProperty(
            material_data_manager.s0,
            "CharacteristicLength",
            245.0
    );

    mgis::behaviour::setMaterialProperty(
            material_data_manager.s1,
            "CharacteristicLength",
            245.0
    );

    mgis::behaviour::setMaterialProperty(
            material_data_manager.s0,
            "PenalisationFactor",
            245.0
    );

    mgis::behaviour::setMaterialProperty(
            material_data_manager.s1,
            "PenalisationFactor",
            245.0
    );

//    mgis::behaviour::setExternalStateVariable(
//            material_data_manager.s0,
//            "FractureEnergy",
//            0.0
//    );
//    mgis::behaviour::setExternalStateVariable(
//            material_data_manager.s1,
//            "FractureEnergy",
//            0.0
//    );

    for (int i = 0; i < material_data_manager.s0.gradients.size(); ++i) {
        print(i, " :" , material_data_manager.s0.gradients[i]);
    }
//    matrix::MatMap<Vector<Real, 5>>(material_data_manager.s0.gradients.data()) = Vector<Real, 5>::Zero();
//    for (int i = 0; i < material_data_manager.s0.gradients.size(); ++i) {
//        print(i, " :" , material_data_manager.s0.gradients[i]);
//    }

    Intg res = mgis::behaviour::integrate(
            material_data_manager,
            mgis::behaviour::IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR,
            0.0,
            0,
            2
    );

    matrix::MatMap<Vector<Real, 5>>(material_data_manager.s0.gradients.data()) = Vector<Real, 5>::Zero();
    for (int i = 0; i < material_data_manager.K.size(); ++i) {
        print(i, " :" , material_data_manager.K[i]);
    }

//    material_data_manager.s0.gradients = Vector<Real, 5>::Zero();
//    material_data_manager.s0.set = Vector<Real, 5>::Zero();
//    print(material_data_manager.s0.gradients);
//    mgis::behaviour::integrate()
//    Intg res = mgis::behaviour::integrate(
//            material_data_manager,
//            mgis::behaviour::IntegrationType::INTEGRATION_CONSISTENT_TANGENT_OPERATOR,
//            0.0,
//            0,
//            1
//    );


}

