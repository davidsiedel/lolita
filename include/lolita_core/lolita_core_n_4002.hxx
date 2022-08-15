#ifndef B9B48BEA_09F5_41DB_84A1_45A6E708901C
#define B9B48BEA_09F5_41DB_84A1_45A6E708901C

#include "lolita_core/lolita.hxx"
#include "lolita_core/lolita_core_n_0000.hxx"
#include "lolita_core/lolita_core_n_1000.hxx"
#include "lolita_core/lolita_core_n_2000.hxx"
#include "lolita_core/lolita_core_n_3000.hxx"
#include "lolita_core/lolita_core_n_4001.hxx"

namespace lolita
{

    template<Element t_element, Domain t_domain>
    struct FiniteElementDegreeOfFreedom
    {

        template<Field t_field, Basis t_basis>
        static constexpr
        Integer
        getSize()
        {
            auto constexpr t_field_size = FieldTraits<t_field>::template getSize<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            return t_field_size * t_basis_size;
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<Real, getSize<t_field, t_basis>()> const>
        getCoefficients()
        const
        {
            auto const & data = degree_of_freedom_->coefficients_.data() + index_;
            return lolita::matrix::Span<lolita::matrix::Vector<Real, getSize<t_field, t_basis>()> const>(data);
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<Real, getSize<t_field, t_basis>()>>
        getCoefficients()
        {
            auto const & data = degree_of_freedom_->coefficients_.data() + index_;
            return lolita::matrix::Span<lolita::matrix::Vector<Real, getSize<t_field, t_basis>()>>(data);
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>
        getCoefficients(
            Integer row,
            Integer col
        )
        {
            auto constexpr t_field_shape = FieldTraits<t_field>::template shape<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            auto const & data = degree_of_freedom_->coefficients_.data() + index_ + (t_field_shape.cols() * row  + col) * t_basis_size;
            return lolita::matrix::Span<lolita::matrix::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>(data);
        }

        template<Field t_field, Basis t_basis>
        lolita::matrix::Span<lolita::matrix::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()> const>
        getCoefficients(
            Integer row,
            Integer col
        )
        const
        {
            auto constexpr t_field_shape = FieldTraits<t_field>::template shape<t_domain>();
            auto constexpr t_basis_size = FiniteElementBasisTraits<t_basis>::template getSize<t_element>();
            auto const & data = degree_of_freedom_->coefficients_.data() + index_ + (t_field_shape.cols() * row  + col) * t_basis_size;
            return lolita::matrix::Span<lolita::matrix::Vector<Real, FiniteElementBasisTraits<t_basis>::template getSize<t_element>()>>(data);
        }
        
        std::shared_ptr<DegreeOfFreedom> const &
        getDegreeOfFreedom()
        const
        {
            return degree_of_freedom_;
        }

        Integer index_;

        std::shared_ptr<DegreeOfFreedom> degree_of_freedom_;

    };

    struct Load
    {

        Load(
            Loading const & loading,
            Integer row,
            Integer col
        )
        :
        row_(row),
        col_(col),
        loading_(loading)
        {}

        Load(
            Loading && loading,
            Integer row,
            Integer col
        )
        :
        row_(row),
        col_(col),
        loading_(std::forward<Loading>(loading))
        {}
        
        Real
        getImposedValue(
            Point const & point,
            Real const & time
        )
        const
        {
            return loading_(point, time);
        }

        Integer
        getRow()
        const
        {
            return row_;
        }

        Integer
        getCol()
        const
        {
            return col_;
        }

        Loading loading_;

        Integer row_;

        Integer col_;

    };
    
}

#endif /* B9B48BEA_09F5_41DB_84A1_45A6E708901C */
