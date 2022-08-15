//
// Created by dsiedel on 12/02/2022.
//

#ifndef FETA_FETA_CORE_BEHAVIOUR_HXX
#define FETA_FETA_CORE_BEHAVIOUR_HXX

#include "new/_feta.hxx"
#include "new/_feta_raise.hxx"
#include "new/_feta_shared_pointer.hxx"

#include <MGIS/Behaviour/Behaviour.hxx>
#include <MGIS/Behaviour/Behaviour.h>
#include <MGIS/Behaviour/MaterialStateManager.h>
#include <MGIS/Behaviour/MaterialStateManager.hxx>
#include <MGIS/Behaviour/Integrate.hxx>
#include <MGIS/Behaviour/Integrate.h>

namespace feta::core
{

    /**
     * @brief
     */
    enum struct HypothesisType
    {

        PlaneStrain,
        PlaneStress,
        AxySymmetric,
        Volumetric

    };

    /**
     * @brief
     */
    enum struct StrainType
    {

        LargeStrain,
        SmallStrain

    };

    /**
     * @brief
     */
    enum struct OutputResult
    {

        BehaviourLawIntegrationFailure,
        BehaviourLawIntegrationSuccess

    };

    /**
     * @brief
     */
    struct Behaviour
    {

        /**
         * @brief
         */
        using BehaviourPointer = SharedPointer<
                mgis::behaviour::Behaviour
        >;

        Behaviour()
        :
                ptr_behaviour(BehaviourPointer()),
                strain_type(StrainType::SmallStrain)
        {}

        Behaviour(
                Strg const &
                path,
                Strg const &
                name,
                HypothesisType const &
                hypothesis_type,
                StrainType const &
                strain_type_arg = StrainType::SmallStrain
        )
        :
                strain_type(strain_type_arg),
                ptr_behaviour(
                getMgisBehaviour(
                        path,
                        name,
                        hypothesis_type,
                        strain_type_arg
                )
        )
        {}

    private:

        /**
         * @brief
         * @param hypothesis_type_arg
         * @return
         */
        inline static
        mgis::behaviour::Hypothesis
        getMgisHypothesis(
                HypothesisType
                hypothesis_type_arg
        )
        {
            mgis::behaviour::Hypothesis hypothesis;
            switch (hypothesis_type_arg) {
                case HypothesisType::PlaneStrain:
                    hypothesis = mgis::behaviour::Hypothesis::PLANESTRAIN;
                    break;
                case HypothesisType::PlaneStress:
                    hypothesis = mgis::behaviour::Hypothesis::PLANESTRESS;
                    break;
                case HypothesisType::AxySymmetric:
                    hypothesis = mgis::behaviour::Hypothesis::AXISYMMETRICAL;
                    break;
                case HypothesisType::Volumetric:
                    hypothesis = mgis::behaviour::Hypothesis::TRIDIMENSIONAL;
                    break;
                default:
                    aassert(false, "unknown HypothesisType");
            }
            return hypothesis;
        }

        static inline
        BehaviourPointer
        getMgisBehaviour(
                Strg const &
                path_arg,
                Strg const &
                name_arg,
                HypothesisType
                hypothesis_type_arg,
                StrainType
                strain_type_arg
        )
        {
            BehaviourPointer a;
            mgis::behaviour::Hypothesis hyp = getMgisHypothesis(
                    hypothesis_type_arg
            );
            if (strain_type_arg == StrainType::SmallStrain) {
                a = BehaviourPointer(
                        mgis::behaviour::load(
                                path_arg,
                                name_arg,
                                hyp
                        )
                );
                return a;
            } else {
                mgis::behaviour::FiniteStrainBehaviourOptions opt;
                opt.stress_measure = mgis::behaviour::FiniteStrainBehaviourOptions::PK1;
                opt.tangent_operator = mgis::behaviour::FiniteStrainBehaviourOptions::TangentOperator::DPK1_DF;
                a = BehaviourPointer(
                        mgis::behaviour::load(
                                opt,
                                path_arg,
                                name_arg,
                                hyp
                        )
                );
                return a;
            }
        }

    public:

        /**
         * @brief
         */
        BehaviourPointer ptr_behaviour;

        /**
         * @brief
         */
        StrainType strain_type;

    };

}

#endif //FETA_FETA_CORE_BEHAVIOUR_HXX
