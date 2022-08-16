#ifndef E11520F9_5899_4E9B_8D95_F11263F73062
#define E11520F9_5899_4E9B_8D95_F11263F73062

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"

namespace lolita
{

    struct DegreeOfFreedom
    {

        DegreeOfFreedom(
            ElementType element_type,
            std::basic_string_view<Character> label
        )
        :
        element_type_(element_type),
        label_(label)
        {}
        
        inline
        Boolean
        operator==(
            DegreeOfFreedom const & other
        )
        const = default;
        
        inline
        Boolean
        operator!=(
            DegreeOfFreedom const & other
        )
        const = default;
        
        inline
        std::basic_string_view<Character>
        getLabel()
        const
        {
            return label_;
        }
        
        inline
        lolita::algebra::Vector<Real> const &
        getCoefficients()
        const
        {
            return coefficients_;
        }
        
        inline
        lolita::algebra::Vector<Real> &
        getCoefficients()
        {
            return coefficients_;
        }

        std::basic_string_view<Character> label_;

        Integer tag_;

        ElementType element_type_;

        lolita::algebra::Vector<Real> coefficients_;

    };

    struct Load2
    {

        Load2(
            ElementType element_type,
            std::basic_string_view<Character> label,
            Loading const & loading,
            Integer row,
            Integer col
        )
        :
        element_type_(element_type),
        label_(label),
        row_(row),
        col_(col),
        loading_(loading)
        {}

        Load2(
            ElementType element_type,
            std::basic_string_view<Character> label,
            Loading && loading,
            Integer row,
            Integer col
        )
        :
        element_type_(element_type),
        label_(label),
        row_(row),
        col_(col),
        loading_(std::forward<Loading>(loading))
        {}
        
        inline
        Boolean
        operator==(
            Load2 const & other
        )
        const = default;
        
        inline
        Boolean
        operator!=(
            Load2 const & other
        )
        const = default;
        
        inline
        std::basic_string_view<Character>
        getLabel()
        const
        {
            return label_;
        }
        
        inline
        Real
        getImposedValue(
            Point const & point,
            Real const & time
        )
        const
        {
            return loading_(point, time);
        }
        
        inline
        Integer
        getRow()
        const
        {
            return row_;
        }
        
        inline
        Integer
        getCol()
        const
        {
            return col_;
        }

        std::basic_string_view<Character> label_;

        ElementType element_type_;

        Integer row_;

        Integer col_;

        Loading loading_;

    };

} // namespace lolita

#endif /* E11520F9_5899_4E9B_8D95_F11263F73062 */
