/**
 * @file region.hxx
 * @author David Siedel (davidsiedel@gmail.com)
 * @brief 
 * @version 0.1
 * @date 2023-03-29
 * 
 * @copyright Copyright (c) 2023
 * 
 */

#ifndef C9EFF827_27A8_4E45_B262_03D5FAF3B724
#define C9EFF827_27A8_4E45_B262_03D5FAF3B724

#include "geometry/frame.hxx"
#include "geometry/shape.hxx"
#include "geometry/domain.hxx"
#include "geometry/point.hxx"

namespace lolita::mesh
{
    
    template<geometry::DomainConcept Domain_, geometry::FrameConcept Frame_>
    struct Region
    {

        explicit
        Region(
            std::basic_string<Character> const & tag
        )
        :
        tag_(tag)
        {}

        explicit
        Region(
            std::basic_string<Character> && tag
        )
        :
        tag_(std::move(tag))
        {}

        std::basic_string<Character> const &
        getLabel()
        const
        {
            return tag_;
        }

    private:

        std::basic_string<Character> tag_;

    };

} // namespace lolita::mesh


#endif /* C9EFF827_27A8_4E45_B262_03D5FAF3B724 */
