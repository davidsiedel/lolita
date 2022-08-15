//
// Created by dsiedel on 12/02/2022.
//

#ifndef FETA_FETA_CORE_MESH_INTERIOR_DOMAIN_HXX
#define FETA_FETA_CORE_MESH_INTERIOR_DOMAIN_HXX

#include "new/_feta_unordered_map.hxx"
#include "new/_feta_shared_pointer.hxx"
//#include "feta/core/element_final.hxx"
#include "new/feta_core_elements.hxx"
#include "new/feta_core_model_description.hxx"
#include "new/feta_core_behaviour.hxx"
#include "new/feta_core_load.hxx"

namespace feta::core
{

//    /**
//     * @brief
//     * @tparam D
//     * @tparam Et
//     */
//    template<Indx D, ElementDescription Et>
//    struct ElementConnectivity;

    namespace element
    {

        template<ElementDescription E, auto M>
        struct Element;

    }

//    /**
//     * @brief
//     * @tparam ElementTypeArg
//     * @tparam FiniteElementDescriptionArg
//     * @tparam MixedElementDescriptionArg
//     */
//    template<ElementDescription ElementTypeArg, auto FiniteElementDescriptionArg, auto MixedElementDescriptionArg>
//    struct FiniteElementCell;

    template<auto Med>
    struct MeshInteriorDomain
    {

        template<elt::ElementDescription Et>
        using Element = elt::Element<Et, Med>;

        using Elements = feta::core::Elements<Med.dim_euclidean, Element>;

    private:

        /**
         * @brief
         */
        using LoadsT = Load3<Med>;

        /**
         * @brief
         */
        using BehaviourPointer = SharedPointer<Behaviour>;

//        /**
//         * @brief
//         */
//        template<ElementDescription Et>
//        using ElementConnectivityPointers = Array<SharedPointer<CellMixedElement<Et, Med>>>;
//
//        /**
//         * @brief
//         */
//        using CellsT = Cells<
//                MixedElementDescriptionArg.getDimEuclidean(),
//                ElementConnectivityPointers
//        >;

    public:

        /**
         * @brief
         */
        MeshInteriorDomain()
        :
        label(""),
        domain_type(elt::Body::Cell),
        loads(""),
        ptr_behaviour()
//        ,
//        cells(CellsT())
        {}

        /**
         * @brief
         * @param label_arg
         * @param path
         * @param name
         * @param hypothesis_type
         * @param strain_type_arg
         * @param load_args
         */
        MeshInteriorDomain(
                Strg &&
                label_arg,
                Strg &&
                path,
                Strg &&
                name,
                HypothesisType &&
                hypothesis_type,
                StrainType &&
                strain_type_arg,
                auto && ...
                load_args
        )
        :
        label(label_arg),
        domain_type(elt::Body::Cell),
        loads(LoadsT{load_args...}),
        ptr_behaviour(Behaviour(path, name, hypothesis_type, strain_type_arg))
//        ,
//        cells(CellsT())
        {}

        /**
         * @brief
         * @param label_arg
         * @param load_args
         */
        MeshInteriorDomain(
                Strg &&
                label_arg,
                elt::Body
                domain_type_arg,
                auto && ...
                load_args
        )
        :
        label(label_arg),
        domain_type(domain_type_arg),
        loads(LoadsT{load_args...}),
        ptr_behaviour()
//        ,
//        cells(CellsT())
        {}

        /**
         * @brief
         */
        elt::Body domain_type;

        /**
         * @brief
         */
        Strg label;

        /**
         * @brief
         */
//        CellsT cells;

        /**
         * @brief
         */
        BehaviourPointer ptr_behaviour;

        /**
         * @brief
         */
        LoadsT loads;

    };

    template<
            auto MixedElementDescriptionArg
    >
    using MeshInteriorDomainPointer = SharedPointer<
            MeshInteriorDomain<
            MixedElementDescriptionArg
    >
    >;

    template<
            auto MixedElementDescriptionArg
    >
    using MeshInteriorDomains = UnorderedMap<
            Strg,
            MeshInteriorDomainPointer<
                    MixedElementDescriptionArg
            >
    >;

//    template<
//            auto FiniteElementDescriptionArg
////            auto MixedElementDescriptionArg
//    >
//    struct MeshInteriorDomainComponent
//    {
//
//    private:
//
//        using LoadT = Load1<
//                FiniteElementDescriptionArg
//        >;
//
////        /**
////         * @brief
////         */
////        template<
////                ElementType Et
////        >
////        using ElementConnectivityPointers = DynamicArray<
////                SharedPointer<
////                        FiniteElementCell<
////                                Et,
////                                FiniteElementDescriptionArg
//////                                MixedElementDescriptionArg
////                        >
////                >
////        >;
////
////        /**
////         * @brief
////         */
////        using Cells = Elements<
////                FiniteElementDescriptionArg.getDimEuclidean(),
////                ElementConnectivityPointers
////        >;
//
//    public:
//
////        MeshInteriorDomain1()
//
//        MeshInteriorDomainComponent()
//        :
//        label(""),
//        load("")
////        cells(Cells())
//        {}
//
//        MeshInteriorDomainComponent(
//                Strg &&
//                label_arg,
//                LoadT const &
//                load_arg,
//                Strg &&
//                path,
//                Strg &&
//                name,
//                HypothesisType &&
//                hypothesis_type,
//                StrainType &&
//                strain_type_arg = StrainType::SmallStrain
//        )
//        :
//        label(label_arg),
//        load(load_arg),
//        behaviour(path, name, hypothesis_type, strain_type_arg)
////        cells(Cells())
//        {}
//
//        /**
//         * @brief
//         */
//        Strg label;
//
//        /**
//         * @brief
//         */
////        Cells cells;
//
//        /**
//         * @brief
//         */
//        Behaviour behaviour;
//
//        /**
//         * @brief
//         */
//        LoadT load;
//
//    };
//
//    template<
//            auto FiniteElementDescriptionArg
////            auto MixedElementDescriptionArg
//    >
//    using MeshInteriorDomainComponentPointer = SharedPointer<
//            MeshInteriorDomainComponent<
//                    FiniteElementDescriptionArg
////                    MixedElementDescriptionArg
//            >
//    >;
//
//
//
////    template<
////            auto MixedElementDescriptionArg
////    >
////    struct ICIQUOI
////    {
////        template<
////                auto FiniteElementDescriptionArg
////        >
////        using MeshInteriorDomainComponentPointer = SharedPointer<
////                MeshInteriorDomainComponent<
////                        FiniteElementDescriptionArg,
////                        MixedElementDescriptionArg
////                >
////        >;
////    };
//
//    template<
//            auto MixedElementDescriptionArg
//    >
//    using MeshInteriorDomain = typename StructuralType<
//            MixedElementDescriptionArg
//    >::template Wrap<
////            MeshInteriorDomain1
//            MeshInteriorDomainComponentPointer
////            typename ICIQUOI<MixedElementDescriptionArg>::MeshInteriorDomainComponentPointer
//    >;

//    template<
//            auto MixedElementDescriptionArg
//    >
//    struct MeshInteriorDomain3
//    :
//    public
//    StructuralType<MixedElementDescriptionArg>::template Wrap<MeshInteriorDomain1>
//    {
//
////        MeshInteriorDomain3(
////                Strg &&
////                label_arg,
////                LoadT const &
////                load_arg,
////                Strg &&
////                path,
////                Strg &&
////                name,
////                HypothesisType &&
////                hypothesis_type,
////                StrainType &&
////                strain_type_arg = StrainType::SmallStrain
////        )
////                :
////                label(label_arg),
////                load(load_arg),
////                behaviour(path, name, hypothesis_type, strain_type_arg)
////        {}
//
//    };

//    /**
//     * @brief
//     * @tparam D
//     */
//    template<
//            Indx D,
//            FieldType F
//    >
//    struct MeshInteriorDomain
//    {
//
//    private:
//
//        using LoadYY = Load<
//                D,
//                F
//        >;
//
//    public:
//
//        /**
//         * @brief
//         */
//        template<
//                ElementType Et
//        >
//        using ElementConnectivityPointers = DynamicArray<
//                SharedPntr<
//                        ElementConnectivity<
//                                D,
//                                Et
//                        >
//                >
//        >;
//
//        /**
//         * @brief
//         */
//        using Cells = Elements<
//                D,
//                ElementConnectivityPointers
//        >;
//
//        /**
//         * @brief
//         */
//        Strg label;
//
//        /**
//         * @brief
//         */
//        Cells cells;
//
//        /**
//         * @brief
//         */
//        Behaviour behaviour;
//
//        /**
//         * @brief
//         */
//        LoadYY load;
//
//    };

}

#endif //FETA_FETA_CORE_MESH_INTERIOR_DOMAIN_HXX
