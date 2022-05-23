//
// Created by dsiedel on 11/05/22.
//

#ifndef LOLITA_LOLITA_CORE_SYSTEM_HXX
#define LOLITA_LOLITA_CORE_SYSTEM_HXX

#include "lolita/lolita.hxx"
#include "lolita/lolita_utility.hxx"
#include "lolita/lolita_algebra.hxx"
#include "lolita/lolita_user.hxx"
#include "lolita/lolita_element.hxx"
#include "lolita/lolita_core_finite_element.hxx"
#include "lolita/lolita_core_mesh.hxx"

namespace lolita::core::mesh
{

//    struct DegreeOfFreedomIndex
//    {
//
//        constexpr
//        DegreeOfFreedomIndex()
//        :
//        unknown_index_(0),
//        binding_index_(0)
//        {}
//
//        constexpr
//        DegreeOfFreedomIndex(
//                lolita::natural &&
//                unknown_index_arg,
//                lolita::natural &&
//                binding_index_arg
//        )
//        :
//        unknown_index_(std::forward<lolita::natural>(unknown_index_arg)),
//        binding_index_(std::forward<lolita::natural>(binding_index_arg))
//        {}
//
//        lolita::natural unknown_index_;
//
//        lolita::natural binding_index_;
//
//    };

    template<lolita::mesh::MeshFormatType Mft, lolita::geometry::Domain _domain, lolita::finite_element::FiniteElementConcept auto... _finite_element>
    struct System
    {

    };

}

#endif //LOLITA_LOLITA_CORE_SYSTEM_HXX
