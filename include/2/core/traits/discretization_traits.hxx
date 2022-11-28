#ifndef C80DDC88_9D89_4048_B55E_126CD1683245
#define C80DDC88_9D89_4048_B55E_126CD1683245

#include "2/core/traits/_include.hxx"
#include "2/core/traits/field_traits.hxx"
#include "2/core/traits/basis_traits.hxx"
#include "2/core/traits/shape_traits.hxx"

namespace lolita::core
{

    template<ElementConcept auto element_>
    struct DiscretizationTraits;

    template<ElementConcept auto element_>
    struct DiscretizationTraitsBase
    {

    private:

        using Implementation_ = DiscretizationTraits<element_>;

    public:

        template<LagrangeShapeConcept auto shape_, MeshConcept auto t_domain>
        static constexpr
        Integer
        getScalarBasisSize()
        {
            return Implementation_::template getScalarBasisSize<shape_, t_domain>();
        }

        template<LagrangeShapeConcept auto shape_, MeshConcept auto t_domain>
        static constexpr
        Integer
        getTensorBasisSize()
        {
            return getScalarBasisSize<shape_, t_domain>() * FieldTraits<element_>::template getSize<t_domain>();
        }

        template<LagrangeShapeConcept auto shape_, MeshConcept auto t_domain>
        static constexpr
        Integer
        getScalarSpaceSize()
        {
            auto num_unknowns = getScalarBasisSize<shape_, t_domain>();
            auto set_num_unknowns = [&] <Integer t_i = 0, Integer t_j = 0> (
                auto & t_set_num_unknowns
            )
            constexpr mutable
            {
                auto constexpr inner_neighbor = ShapeTraits<shape_>::template getInnerNeighbor<t_i, t_j>();
                auto constexpr num_inner_neighbors = ShapeTraits<shape_>::template getNumInnerNeighbors<t_i, t_j>();
                num_unknowns += getScalarBasisSize<inner_neighbor, t_domain>() * num_inner_neighbors;
                if constexpr (t_j < ShapeTraits<shape_>::template getNumInnerNeighbors<t_i>() - 1)
                {
                    t_set_num_unknowns.template operator ()<t_i, t_j + 1>(t_set_num_unknowns);
                }
                else if constexpr (t_i < ShapeTraits<shape_>::template getNumInnerNeighbors<>() - 1)
                {
                    t_set_num_unknowns.template operator ()<t_i + 1, 0>(t_set_num_unknowns);
                }
            };
            if constexpr (shape_ != Node())
            {
                set_num_unknowns(set_num_unknowns);
            }
            return num_unknowns;
        }

        template<LagrangeShapeConcept auto shape_, MeshConcept auto t_domain>
        static constexpr
        Integer
        getTensorSpaceSize()
        {
            return getScalarSpaceSize<shape_, t_domain>() * FieldTraits<element_>::template getSize<t_domain>();
        }

    };

    template<HdgElementConcept auto element_>
    struct DiscretizationTraits<element_> : DiscretizationTraitsBase<element_>
    {

        template<LagrangeShapeConcept auto shape_, MeshConcept auto mesh_>
        static constexpr
        Integer
        getScalarBasisSize()
        {
            return 0;
        }

        template<LagrangeShapeConcept auto shape_, MeshConcept auto mesh_>
        static constexpr
        Integer
        getScalarBasisSize()
        requires(shape_.getDim() == element_.getDimDomain())
        {
            return BasisTraits<element_.getCellBasis()>::template getSize<shape_>();
        }

        template<LagrangeShapeConcept auto shape_, MeshConcept auto mesh_>
        static constexpr
        Integer
        getScalarBasisSize()
        requires(shape_.getDim() == element_.getDimDomain() - 1)
        {
            return BasisTraits<element_.getFaceBasis()>::template getSize<shape_>();
        }
        
    };

} // namespace lolita::core

#endif /* C80DDC88_9D89_4048_B55E_126CD1683245 */
