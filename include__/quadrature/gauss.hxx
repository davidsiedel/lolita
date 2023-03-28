#ifndef E603E856_5AEE_4096_91AA_872A653000CB
#define E603E856_5AEE_4096_91AA_872A653000CB

#include "geometry/shape.hxx"

namespace lolita::quadrature
{

    template<geometry::ShapeConcept Shape_, Integer k_>
    struct GaussQuadrature;

    template<>
    struct GaussQuadrature<geometry::Node, 0>
    {

        static constexpr
        Integer size_ = 1;
        
        static constexpr 
        std::array<std::array<Real, 3>, size_> reference_points_ = {
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
        };
        
        static constexpr 
        std::array<Real, size_> reference_weights_ = {
            +1.0000000000000000,
        };

    };

    template<>
    struct GaussQuadrature<geometry::Segment, 0>
    {

        static constexpr
        Integer size_ = 1;
        
        static constexpr 
        std::array<std::array<Real, 3>, size_> reference_points_ = {
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
        };
        
        static constexpr 
        std::array<Real, size_> reference_weights_ = {
            +2.0000000000000000,
        };

    };

    template<>
    struct GaussQuadrature<geometry::Triangle, 0>
    {

        static constexpr
        Integer size_ = 1;
        
        static constexpr 
        std::array<std::array<Real, 3>, size_> reference_points_ = {
            +0.3333333333333333, +0.3333333333333333, +0.0000000000000000,
        };
        
        static constexpr 
        std::array<Real, size_> reference_weights_ = {
            +0.5000000000000000,
        };

    };
    
} // namespace lolita::quadrature


#endif /* E603E856_5AEE_4096_91AA_872A653000CB */
