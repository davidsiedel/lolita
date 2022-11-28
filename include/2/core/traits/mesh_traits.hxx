#ifndef CCBCD4E5_FF13_40E7_BC77_E209E3883914
#define CCBCD4E5_FF13_40E7_BC77_E209E3883914

#include "2/core/traits/_include.hxx"
#include "2/core/traits/shape_traits.hxx"
#include "2/core/traits/domain_traits.hxx"

namespace lolita::core
{
    
    template<MeshConcept auto t_domain>
    struct MeshTraits
    {

    private:

        using Elements_ = lolita::utility::tuple_slice_t<typename ShapeLibrary::Elements<ShapeView>, 0, t_domain.getDim() + 1>;

        using Domains_ = lolita::utility::tuple_slice_t<typename DomainLibrary::Domains<DomainView>, 0, t_domain.getDim() + 1>;

    public:

        static constexpr
        Boolean
        hasCartesianCoordinates()
        {
            return CartesianMeshConcept<decltype(t_domain)>;
        }

        static constexpr
        Boolean
        hasAxiSymmetricCoordinates()
        {
            return AxiSymmetricMeshConcept<decltype(t_domain)>;
        }

        static constexpr
        Integer
        letQuadratureOrder(
            Integer polynomial_order
        )
        {
            return hasAxiSymmetricCoordinates() ? 2 * polynomial_order + 1 : 2 * polynomial_order;
        }
        
        template<LagrangeShapeConcept auto t_element>
        static constexpr
        ShapeCoordinates
        getElementCoordinates()
        {
            auto coordinates = ShapeCoordinates(t_element.getDim(), -1);
            auto set_coordinates = [&] <Integer t_i = 0> (
                auto & t_set_coordinates
            )
            constexpr mutable
            {
                if (std::tuple_element_t<t_i, std::tuple_element_t<t_element.getDim(), Elements_>>::getElement() == t_element)
                {
                    coordinates.setTag(t_i);
                }
                if constexpr (t_i < std::tuple_size_v<std::tuple_element_t<t_element.getDim(), Elements_>> - 1)
                {
                    t_set_coordinates.template operator()<t_i + 1>(t_set_coordinates);
                }
            };
            set_coordinates(set_coordinates);
            return coordinates;
        }
        
        template<Integer t_i, Integer t_j>
        static constexpr
        LagrangeShapeConcept auto
        getElement()
        {
            return std::tuple_element_t<t_j, std::tuple_element_t<t_i, Elements_>>::getElement();
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumElements()
        requires(sizeof...(t_i) == 0)
        {
            return std::tuple_size_v<Elements_>;
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumElements()
        requires(sizeof...(t_i) == 1)
        {
            return std::tuple_size_v<std::tuple_element_t<std::array<Integer, 1>{t_i...}[0], Elements_>>;
        }
        
        template<Integer t_i>
        static constexpr
        DomainConcept auto const &
        getDomain()
        {
            return std::tuple_element_t<t_i, Domains_>::getDomain();
        }
        
        template<Integer... t_i>
        static constexpr
        Integer
        getNumDomains()
        requires(sizeof...(t_i) == 0)
        {
            return std::tuple_size_v<Domains_>;
        }

    };
    
} // namespace lolita::core


#endif /* CCBCD4E5_FF13_40E7_BC77_E209E3883914 */
