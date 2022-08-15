//
// Created by dsiedel on 09/02/2022.
//

#ifndef FETA_FETA_CORE_ELEMENTS_HXX
#define FETA_FETA_CORE_ELEMENTS_HXX

#include "new/_feta_aggregate.hxx"
#include "new/_feta_collection.hxx"
#include "new/feta_core_element_element_description.hxx"

namespace feta::core
{

    namespace elt = feta::core::element;

    template<template<elt::ElementDescription> typename T>
    using PointElements = Collection<T<elt::pnt_0>>;

    template<template<elt::ElementDescription> typename T>
    using CurveElements = Collection<T<elt::seg_2>>;

    template<template<elt::ElementDescription> typename T>
    using PlaneElements = Collection<T<elt::tri_3>, T<elt::qua_4>>;

    template<template<elt::ElementDescription> typename T>
    using SolidElements = Collection<T<elt::tet_4>>;

    namespace detail
    {

        template<Indx D, template<elt::ElementDescription> typename T>
        struct ElementsPolicy;

        template<template<elt::ElementDescription> typename T>
        struct ElementsPolicy<1, T>
        {

            using Type = Collection<PointElements<T>, CurveElements<T>>;

        };

        template<template<elt::ElementDescription> typename T>
        struct ElementsPolicy<2, T>
        {

            using Type = Collection<PointElements<T>, CurveElements<T>, PlaneElements<T>>;

        };

        template<template<elt::ElementDescription> typename T>
        struct ElementsPolicy<3, T>
        {

            using Type = Collection<PointElements<T>, CurveElements<T>, PlaneElements<T>, SolidElements<T>>;

        };

    }

    /**
     * @brief
     */
    template<Indx D, template<elt::ElementDescription> typename T>
    using Elements = typename detail::ElementsPolicy<D, T>::Type;

}

#endif //FETA_FETA_CORE_ELEMENTS_HXX
