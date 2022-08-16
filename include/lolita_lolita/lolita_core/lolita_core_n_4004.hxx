#ifndef CB51704E_040A_4B3D_AE3C_46C0DE3B8543
#define CB51704E_040A_4B3D_AE3C_46C0DE3B8543

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_3000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4001.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4002.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_4003.hxx"

namespace lolita
{

    template<Element t_element, Domain t_domain>
    struct FiniteElement
    {

        std::shared_ptr<std::vector<QuadraturePoint2>> element_integration_points2_;

        std::shared_ptr<ElementIntegrationPoints<t_domain>> element_integration_points_;

        // std::vector<std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>>> degrees_of_freedom_;

        std::vector<std::shared_ptr<mgis::behaviour::BehaviourData>> behavior_data_;

        std::vector<lolita::algebra::Matrix<Real>> element_operators_;

        // std::vector<std::shared_ptr<Load>> loads_;
        
        std::shared_ptr<mgis::behaviour::Behaviour> behavior_;

        std::shared_ptr<std::basic_string<Character>> label_;

        FiniteElement(
            std::shared_ptr<std::basic_string<Character>> const & label
        )
        :
        label_(label)
        {}

        std::basic_string_view<Character>
        getLabel()
        const
        {
            return std::basic_string_view<Character>(* label_);
        }

        // Boolean
        // hasLoad(
        //     Integer row,
        //     Integer col
        // )
        // const
        // {
        //     auto has_load = [&] (
        //         std::shared_ptr<Load> const & load
        //     )
        //     {
        //         return load->getRow() == row && load->getCol() == col;
        //     };
        //     return std::find_if(loads_.begin(), loads_.end(), has_load) != loads_.end();
        // }

        // Integer
        // getLoadIndex(
        //     Integer row,
        //     Integer col
        // )
        // const
        // {
        //     auto has_load = [&] (
        //         std::shared_ptr<Load> const & load
        //     )
        //     {
        //         return load->getRow() == row && load->getCol() == col;
        //     };
        //     if (std::find_if(loads_.begin(), loads_.end(), has_load) != loads_.end())
        //     {
        //         return std::distance(loads_.begin(), std::find_if(loads_.begin(), loads_.end(), has_load));
        //     }
        //     else
        //     {
        //         throw std::runtime_error("NO");
        //     }
        // }
        
        // void
        // addLoad(
        //     std::shared_ptr<Load> const & load
        // )
        // {
        //     if (!hasLoad(load->getRow(), load->getCol()))
        //     {
        //         loads_.push_back(load);
        //     }
        //     else
        //     {
        //         getLoad(load->getRow(), load->getCol()) = load;
        //     }
        // }

        // Boolean
        // hasDegreeOfFreedom(
        //     std::basic_string_view<Character> label
        // )
        // const
        // {
        //     auto has_label = [&] (
        //         std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const & degree_of_freedom
        //     )
        //     {
        //         return degree_of_freedom->getDegreeOfFreedom()->getLabel() == label;
        //     };
        //     return std::find_if(degrees_of_freedom_.begin(), degrees_of_freedom_.end(), has_label) != degrees_of_freedom_.end();
        // }

        // Integer
        // getDegreeOfFreedomIndex(
        //     std::basic_string_view<Character> label
        // )
        // const
        // {
        //     auto has_label = [&] (
        //         std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const & degree_of_freedom
        //     )
        //     {
        //         return degree_of_freedom->getDegreeOfFreedom()->getLabel() == label;
        //     };
        //     if (std::find_if(degrees_of_freedom_.begin(), degrees_of_freedom_.end(), has_label) != degrees_of_freedom_.end())
        //     {
        //         return std::distance(degrees_of_freedom_.begin(), std::find_if(degrees_of_freedom_.begin(), degrees_of_freedom_.end(), has_label));
        //     }
        //     else
        //     {
        //         throw std::runtime_error("NO");
        //     }
        // }

        // template<Field t_field, Basis t_basis>
        // void
        // addDegreeOfFreedom(
        //     std::shared_ptr<DegreeOfFreedom> & degree_of_freedom
        // )
        // {
        //     auto constexpr t_size = FiniteElementDegreeOfFreedom<t_element, t_domain>::template getSize<t_field, t_basis>();
        //     if (!hasDegreeOfFreedom(degree_of_freedom->getLabel()))
        //     {
        //         auto element_dof = std::make_unique<FiniteElementDegreeOfFreedom<t_element, t_domain>>();
        //         element_dof->index_ = degree_of_freedom->coefficients_.size();
        //         element_dof->degree_of_freedom_ = degree_of_freedom;
        //         degree_of_freedom->coefficients_.resize(degree_of_freedom->coefficients_.size() + t_size);
        //         degrees_of_freedom_.push_back(std::move(element_dof));
        //     }
        // }

        // std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> const &
        // getDegreeOfFreedom(
        //     std::basic_string_view<Character> label
        // )
        // const
        // {
        //     return degrees_of_freedom_[getDegreeOfFreedomIndex(label)];
        // }

        // std::unique_ptr<FiniteElementDegreeOfFreedom<t_element, t_domain>> &
        // getDegreeOfFreedom(
        //     std::basic_string_view<Character> label
        // )
        // {
        //     return degrees_of_freedom_[getDegreeOfFreedomIndex(label)];
        // }

        // std::shared_ptr<Load> const &
        // getLoad(
        //     Integer row,
        //     Integer col
        // )
        // const
        // {
        //     return loads_[getLoadIndex(row, col)];
        // }

        // std::shared_ptr<Load> &
        // getLoad(
        //     Integer row,
        //     Integer col
        // )
        // {
        //     return loads_[getLoadIndex(row, col)];
        // }

    };
    
} // namespace lolita



#endif /* CB51704E_040A_4B3D_AE3C_46C0DE3B8543 */
