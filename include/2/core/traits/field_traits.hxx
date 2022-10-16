#ifndef A1E393BC_29F1_45CE_8EE0_DD2097797B01
#define A1E393BC_29F1_45CE_8EE0_DD2097797B01

#include "2/core/traits/_include.hxx"

namespace lolita::core
{

    template<FieldConcept auto field_>
    struct FieldView
    {

        static constexpr
        FieldConcept auto const &
        getField()
        {
            return field_;
        }

    };

    template<FieldConcept auto field_>
    struct FieldTraits
    {

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumCols()
        requires(field_.isTensor(0))
        {
            return lolita::numerics::pow(mesh_.getDim(), 0);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumCols()
        requires(field_.isTensor(1))
        {
            return lolita::numerics::pow(mesh_.getDim(), 0);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumCols()
        requires(field_.isTensor(2))
        {
            return lolita::numerics::pow(mesh_.getDim(), 1);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumCols()
        requires(field_.isTensor(3))
        {
            return lolita::numerics::pow(mesh_.getDim(), 1);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumCols()
        requires(field_.isTensor(4))
        {
            return lolita::numerics::pow(mesh_.getDim(), 2);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumRows()
        requires(field_.isTensor(0))
        {
            return lolita::numerics::pow(mesh_.getDim(), 0);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumRows()
        requires(field_.isTensor(1))
        {
            return lolita::numerics::pow(mesh_.getDim(), 1);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumRows()
        requires(field_.isTensor(2))
        {
            return lolita::numerics::pow(mesh_.getDim(), 1);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumRows()
        requires(field_.isTensor(3))
        {
            return lolita::numerics::pow(mesh_.getDim(), 2);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getNumRows()
        requires(field_.isTensor(4))
        {
            return lolita::numerics::pow(mesh_.getDim(), 2);
        }

        template<MeshConcept auto mesh_>
        static constexpr
        Integer
        getSize()
        {
            return getNumRows<mesh_>() * getNumCols<mesh_>();
        }

    };
    
} // namespace lolita::core


#endif /* A1E393BC_29F1_45CE_8EE0_DD2097797B01 */
