#ifndef CC79BDC5_49DB_4A81_8A93_18ABD6551AF1
#define CC79BDC5_49DB_4A81_8A93_18ABD6551AF1

#include "lolita_lolita/lolita_core/lolita.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_0000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_1000.hxx"
#include "lolita_lolita/lolita_core/lolita_core_n_2000.hxx"

namespace lolita
{
    
    template<Element t_element, Quadrature t_quadrature>
    struct ElementQuadratureRuleTraits;
    
    template<Element t_element, Quadrature t_quadrature>
    // requires(t_element.isNode() || !t_element.isNode())
    struct ElementQuadratureRuleTraits
    {

        static constexpr
        Integer
        getSize()
        {
            return 1;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
            +2.0000000000000000
        };

    };
    
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // SEGMENT
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Segment) && t_quadrature.hasOrd(1))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 1;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
            +0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
            +2.0000000000000000,
        };

    };
    
    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Segment) && t_quadrature.hasOrd(2))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 2;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.5773502691896257, +0.0000000000000000, +0.0000000000000000,
+0.5773502691896257, +0.0000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+1.0000000000000000,
+1.0000000000000000,
        };

    };
    
    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Segment) && t_quadrature.hasOrd(3))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 3;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.7745966692414834, +0.0000000000000000, +0.0000000000000000,
+0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
+0.7745966692414834, +0.0000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.5555555555555555,
+0.8888888888888892,
+0.5555555555555555,
        };

    };
    
    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Segment) && t_quadrature.hasOrd(4))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 4;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.8611363115940526, +0.0000000000000000, +0.0000000000000000,
-0.3399810435848563, +0.0000000000000000, +0.0000000000000000,
+0.3399810435848563, +0.0000000000000000, +0.0000000000000000,
+0.8611363115940526, +0.0000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.3478548451374536,
+0.6521451548625466,
+0.6521451548625466,
+0.3478548451374536,
        };

    };
    
    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Segment) && t_quadrature.hasOrd(5))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 5;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.9061798459386641, +0.0000000000000000, +0.0000000000000000,
-0.5384693101056831, +0.0000000000000000, +0.0000000000000000,
+0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
+0.5384693101056831, +0.0000000000000000, +0.0000000000000000,
+0.9061798459386641, +0.0000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.2369268850561891,
+0.4786286704993662,
+0.5688888888888890,
+0.4786286704993662,
+0.2369268850561891,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Segment) && t_quadrature.hasOrd(6))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 6;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.9324695142031521, +0.0000000000000000, +0.0000000000000000,
-0.6612093864662645, +0.0000000000000000, +0.0000000000000000,
-0.2386191860831969, +0.0000000000000000, +0.0000000000000000,
+0.2386191860831969, +0.0000000000000000, +0.0000000000000000,
+0.6612093864662645, +0.0000000000000000, +0.0000000000000000,
+0.9324695142031521, +0.0000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.1713244923791704,
+0.3607615730481387,
+0.4679139345726911,
+0.4679139345726911,
+0.3607615730481387,
+0.1713244923791704,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Segment) && t_quadrature.hasOrd(7))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 7;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.9491079123427585, +0.0000000000000000, +0.0000000000000000,
-0.7415311855993945, +0.0000000000000000, +0.0000000000000000,
-0.4058451513773972, +0.0000000000000000, +0.0000000000000000,
+0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
+0.4058451513773972, +0.0000000000000000, +0.0000000000000000,
+0.7415311855993945, +0.0000000000000000, +0.0000000000000000,
+0.9491079123427585, +0.0000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.1294849661688702,
+0.2797053914892766,
+0.3818300505051188,
+0.4179591836734692,
+0.3818300505051188,
+0.2797053914892766,
+0.1294849661688702,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Segment) && t_quadrature.hasOrd(8))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 8;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.9602898564975362, +0.0000000000000000, +0.0000000000000000,
-0.7966664774136267, +0.0000000000000000, +0.0000000000000000,
-0.5255324099163290, +0.0000000000000000, +0.0000000000000000,
-0.1834346424956498, +0.0000000000000000, +0.0000000000000000,
+0.1834346424956498, +0.0000000000000000, +0.0000000000000000,
+0.5255324099163290, +0.0000000000000000, +0.0000000000000000,
+0.7966664774136267, +0.0000000000000000, +0.0000000000000000,
+0.9602898564975362, +0.0000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.1012285362903764,
+0.2223810344533742,
+0.3137066458778874,
+0.3626837833783620,
+0.3626837833783620,
+0.3137066458778874,
+0.2223810344533742,
+0.1012285362903764,
        };

    };

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // TRIANGLE
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Triangle) && t_quadrature.hasOrd(1))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 1;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
+0.3333333333333333, +0.3333333333333333, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.5000000000000000,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Triangle) && t_quadrature.hasOrd(2))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 3;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
+0.1666666666666667, +0.6666666666666666, +0.0000000000000000,
+0.6666666666666666, +0.1666666666666667, +0.0000000000000000,
+0.1666666666666667, +0.1666666666666667, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.1666666666666667,
+0.1666666666666667,
+0.1666666666666667,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Triangle) && t_quadrature.hasOrd(3))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 4;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
+0.3333333333333333, +0.3333333333333333, +0.0000000000000000,
+0.2000000000000000, +0.6000000000000000, +0.0000000000000000,
+0.6000000000000000, +0.2000000000000000, +0.0000000000000000,
+0.2000000000000000, +0.2000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
-0.2812500000000000,
+0.2604166666666667,
+0.2604166666666667,
+0.2604166666666667,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Triangle) && t_quadrature.hasOrd(4))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 6;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
+0.4459484909159649, +0.1081030181680702, +0.0000000000000000,
+0.0915762135097707, +0.8168475729804585, +0.0000000000000000,
+0.1081030181680702, +0.4459484909159649, +0.0000000000000000,
+0.8168475729804585, +0.0915762135097707, +0.0000000000000000,
+0.4459484909159649, +0.4459484909159649, +0.0000000000000000,
+0.0915762135097707, +0.0915762135097707, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.1116907948390058,
+0.0549758718276609,
+0.1116907948390058,
+0.0549758718276609,
+0.1116907948390058,
+0.0549758718276609,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Triangle) && t_quadrature.hasOrd(5))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 7;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
+0.3333333333333333, +0.3333333333333333, +0.0000000000000000,
+0.4701420641051151, +0.0597158717897699, +0.0000000000000000,
+0.1012865073234563, +0.7974269853530873, +0.0000000000000000,
+0.0597158717897699, +0.4701420641051151, +0.0000000000000000,
+0.7974269853530873, +0.1012865073234563, +0.0000000000000000,
+0.4701420641051151, +0.4701420641051151, +0.0000000000000000,
+0.1012865073234563, +0.1012865073234563, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.1125000000000000,
+0.0661970763942531,
+0.0629695902724136,
+0.0661970763942531,
+0.0629695902724136,
+0.0661970763942531,
+0.0629695902724136,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Triangle) && t_quadrature.hasOrd(6))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 12;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
+0.2492867451709104, +0.5014265096581791, +0.0000000000000000,
+0.0630890144915023, +0.8738219710169954, +0.0000000000000000,
+0.5014265096581791, +0.2492867451709104, +0.0000000000000000,
+0.8738219710169954, +0.0630890144915023, +0.0000000000000000,
+0.2492867451709104, +0.2492867451709104, +0.0000000000000000,
+0.0630890144915023, +0.0630890144915023, +0.0000000000000000,
+0.3103524510337844, +0.6365024991213987, +0.0000000000000000,
+0.0531450498448170, +0.3103524510337844, +0.0000000000000000,
+0.6365024991213987, +0.0531450498448170, +0.0000000000000000,
+0.0531450498448170, +0.6365024991213987, +0.0000000000000000,
+0.3103524510337844, +0.0531450498448170, +0.0000000000000000,
+0.6365024991213987, +0.3103524510337844, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.0583931378631897,
+0.0254224531851035,
+0.0583931378631897,
+0.0254224531851035,
+0.0583931378631897,
+0.0254224531851035,
+0.0414255378091868,
+0.0414255378091868,
+0.0414255378091868,
+0.0414255378091868,
+0.0414255378091868,
+0.0414255378091868,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Triangle) && t_quadrature.hasOrd(7))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 13;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
+0.3333333333333333, +0.3333333333333333, +0.0000000000000000,
+0.2603459660790399, +0.4793080678419203, +0.0000000000000000,
+0.0651301029022158, +0.8697397941955683, +0.0000000000000000,
+0.4793080678419203, +0.2603459660790399, +0.0000000000000000,
+0.8697397941955683, +0.0651301029022158, +0.0000000000000000,
+0.2603459660790399, +0.2603459660790399, +0.0000000000000000,
+0.0651301029022158, +0.0651301029022158, +0.0000000000000000,
+0.3128654960048739, +0.6384441885698097, +0.0000000000000000,
+0.0486903154253164, +0.3128654960048739, +0.0000000000000000,
+0.6384441885698097, +0.0486903154253164, +0.0000000000000000,
+0.0486903154253164, +0.6384441885698097, +0.0000000000000000,
+0.3128654960048739, +0.0486903154253164, +0.0000000000000000,
+0.6384441885698097, +0.3128654960048739, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
-0.0747850222338411,
+0.0878076287166039,
+0.0266736178044192,
+0.0878076287166039,
+0.0266736178044192,
+0.0878076287166039,
+0.0266736178044192,
+0.0385568804451286,
+0.0385568804451286,
+0.0385568804451286,
+0.0385568804451286,
+0.0385568804451286,
+0.0385568804451286,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Triangle) && t_quadrature.hasOrd(8))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 16;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
+0.3333333333333333, +0.3333333333333333, +0.0000000000000000,
+0.4592925882927232, +0.0814148234145536, +0.0000000000000000,
+0.1705693077517603, +0.6588613844964795, +0.0000000000000000,
+0.0505472283170311, +0.8989055433659379, +0.0000000000000000,
+0.0814148234145536, +0.4592925882927232, +0.0000000000000000,
+0.6588613844964795, +0.1705693077517603, +0.0000000000000000,
+0.8989055433659379, +0.0505472283170311, +0.0000000000000000,
+0.4592925882927232, +0.4592925882927232, +0.0000000000000000,
+0.1705693077517603, +0.1705693077517603, +0.0000000000000000,
+0.0505472283170311, +0.0505472283170311, +0.0000000000000000,
+0.2631128296346381, +0.7284923929554042, +0.0000000000000000,
+0.0083947774099577, +0.2631128296346381, +0.0000000000000000,
+0.7284923929554042, +0.0083947774099577, +0.0000000000000000,
+0.0083947774099577, +0.7284923929554042, +0.0000000000000000,
+0.2631128296346381, +0.0083947774099577, +0.0000000000000000,
+0.7284923929554042, +0.2631128296346381, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.0721578038388936,
+0.0475458171336423,
+0.0516086852673591,
+0.0162292488115990,
+0.0475458171336423,
+0.0516086852673591,
+0.0162292488115990,
+0.0475458171336423,
+0.0516086852673591,
+0.0162292488115990,
+0.0136151570872175,
+0.0136151570872175,
+0.0136151570872175,
+0.0136151570872175,
+0.0136151570872175,
+0.0136151570872175,
        };

    };

    // ---------------------------------------------------------------------------------------------------------------------------------------------------------
    // QUADRANGLE
    // ---------------------------------------------------------------------------------------------------------------------------------------------------------

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Quadrangle) && t_quadrature.hasOrd(1))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 1;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
+0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+4.0000000000000000,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Quadrangle) && t_quadrature.hasOrd(2))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 4;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.5773502691896257, -0.5773502691896257, +0.0000000000000000,
+0.5773502691896257, -0.5773502691896257, +0.0000000000000000,
-0.5773502691896257, +0.5773502691896257, +0.0000000000000000,
+0.5773502691896257, +0.5773502691896257, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+1.0000000000000000,
+1.0000000000000000,
+1.0000000000000000,
+1.0000000000000000,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Quadrangle) && t_quadrature.hasOrd(3))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 9;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.7745966692414834, -0.7745966692414834, +0.0000000000000000,
+0.0000000000000000, -0.7745966692414834, +0.0000000000000000,
+0.7745966692414834, -0.7745966692414834, +0.0000000000000000,
-0.7745966692414834, +0.0000000000000000, +0.0000000000000000,
+0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
+0.7745966692414834, +0.0000000000000000, +0.0000000000000000,
-0.7745966692414834, +0.7745966692414834, +0.0000000000000000,
+0.0000000000000000, +0.7745966692414834, +0.0000000000000000,
+0.7745966692414834, +0.7745966692414834, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.3086419753086419,
+0.4938271604938272,
+0.3086419753086419,
+0.4938271604938272,
+0.7901234567901240,
+0.4938271604938272,
+0.3086419753086419,
+0.4938271604938272,
+0.3086419753086419,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Quadrangle) && t_quadrature.hasOrd(4))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 16;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.8611363115940526, -0.8611363115940526, +0.0000000000000000,
-0.3399810435848563, -0.8611363115940526, +0.0000000000000000,
+0.3399810435848563, -0.8611363115940526, +0.0000000000000000,
+0.8611363115940526, -0.8611363115940526, +0.0000000000000000,
-0.8611363115940526, -0.3399810435848563, +0.0000000000000000,
-0.3399810435848563, -0.3399810435848563, +0.0000000000000000,
+0.3399810435848563, -0.3399810435848563, +0.0000000000000000,
+0.8611363115940526, -0.3399810435848563, +0.0000000000000000,
-0.8611363115940526, +0.3399810435848563, +0.0000000000000000,
-0.3399810435848563, +0.3399810435848563, +0.0000000000000000,
+0.3399810435848563, +0.3399810435848563, +0.0000000000000000,
+0.8611363115940526, +0.3399810435848563, +0.0000000000000000,
-0.8611363115940526, +0.8611363115940526, +0.0000000000000000,
-0.3399810435848563, +0.8611363115940526, +0.0000000000000000,
+0.3399810435848563, +0.8611363115940526, +0.0000000000000000,
+0.8611363115940526, +0.8611363115940526, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.1210029932856018,
+0.2268518518518518,
+0.2268518518518518,
+0.1210029932856018,
+0.2268518518518518,
+0.4252933030106950,
+0.4252933030106950,
+0.2268518518518518,
+0.2268518518518518,
+0.4252933030106950,
+0.4252933030106950,
+0.2268518518518518,
+0.1210029932856018,
+0.2268518518518518,
+0.2268518518518518,
+0.1210029932856018,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Quadrangle) && t_quadrature.hasOrd(5))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 25;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.9061798459386641, -0.9061798459386641, +0.0000000000000000,
-0.5384693101056831, -0.9061798459386641, +0.0000000000000000,
+0.0000000000000000, -0.9061798459386641, +0.0000000000000000,
+0.5384693101056831, -0.9061798459386641, +0.0000000000000000,
+0.9061798459386641, -0.9061798459386641, +0.0000000000000000,
-0.9061798459386641, -0.5384693101056831, +0.0000000000000000,
-0.5384693101056831, -0.5384693101056831, +0.0000000000000000,
+0.0000000000000000, -0.5384693101056831, +0.0000000000000000,
+0.5384693101056831, -0.5384693101056831, +0.0000000000000000,
+0.9061798459386641, -0.5384693101056831, +0.0000000000000000,
-0.9061798459386641, +0.0000000000000000, +0.0000000000000000,
-0.5384693101056831, +0.0000000000000000, +0.0000000000000000,
+0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
+0.5384693101056831, +0.0000000000000000, +0.0000000000000000,
+0.9061798459386641, +0.0000000000000000, +0.0000000000000000,
-0.9061798459386641, +0.5384693101056831, +0.0000000000000000,
-0.5384693101056831, +0.5384693101056831, +0.0000000000000000,
+0.0000000000000000, +0.5384693101056831, +0.0000000000000000,
+0.5384693101056831, +0.5384693101056831, +0.0000000000000000,
+0.9061798459386641, +0.5384693101056831, +0.0000000000000000,
-0.9061798459386641, +0.9061798459386641, +0.0000000000000000,
-0.5384693101056831, +0.9061798459386641, +0.0000000000000000,
+0.0000000000000000, +0.9061798459386641, +0.0000000000000000,
+0.5384693101056831, +0.9061798459386641, +0.0000000000000000,
+0.9061798459386641, +0.9061798459386641, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.0561343488624286,
+0.1133999999999999,
+0.1347850723875209,
+0.1133999999999999,
+0.0561343488624286,
+0.1133999999999999,
+0.2290854042239909,
+0.2722865325507506,
+0.2290854042239909,
+0.1133999999999999,
+0.1347850723875209,
+0.2722865325507506,
+0.3236345679012347,
+0.2722865325507506,
+0.1347850723875209,
+0.1133999999999999,
+0.2290854042239909,
+0.2722865325507506,
+0.2290854042239909,
+0.1133999999999999,
+0.0561343488624286,
+0.1133999999999999,
+0.1347850723875209,
+0.1133999999999999,
+0.0561343488624286,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Quadrangle) && t_quadrature.hasOrd(6))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 36;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.9324695142031521, -0.9324695142031521, +0.0000000000000000,
-0.6612093864662645, -0.9324695142031521, +0.0000000000000000,
-0.2386191860831969, -0.9324695142031521, +0.0000000000000000,
+0.2386191860831969, -0.9324695142031521, +0.0000000000000000,
+0.6612093864662645, -0.9324695142031521, +0.0000000000000000,
+0.9324695142031521, -0.9324695142031521, +0.0000000000000000,
-0.9324695142031521, -0.6612093864662645, +0.0000000000000000,
-0.6612093864662645, -0.6612093864662645, +0.0000000000000000,
-0.2386191860831969, -0.6612093864662645, +0.0000000000000000,
+0.2386191860831969, -0.6612093864662645, +0.0000000000000000,
+0.6612093864662645, -0.6612093864662645, +0.0000000000000000,
+0.9324695142031521, -0.6612093864662645, +0.0000000000000000,
-0.9324695142031521, -0.2386191860831969, +0.0000000000000000,
-0.6612093864662645, -0.2386191860831969, +0.0000000000000000,
-0.2386191860831969, -0.2386191860831969, +0.0000000000000000,
+0.2386191860831969, -0.2386191860831969, +0.0000000000000000,
+0.6612093864662645, -0.2386191860831969, +0.0000000000000000,
+0.9324695142031521, -0.2386191860831969, +0.0000000000000000,
-0.9324695142031521, +0.2386191860831969, +0.0000000000000000,
-0.6612093864662645, +0.2386191860831969, +0.0000000000000000,
-0.2386191860831969, +0.2386191860831969, +0.0000000000000000,
+0.2386191860831969, +0.2386191860831969, +0.0000000000000000,
+0.6612093864662645, +0.2386191860831969, +0.0000000000000000,
+0.9324695142031521, +0.2386191860831969, +0.0000000000000000,
-0.9324695142031521, +0.6612093864662645, +0.0000000000000000,
-0.6612093864662645, +0.6612093864662645, +0.0000000000000000,
-0.2386191860831969, +0.6612093864662645, +0.0000000000000000,
+0.2386191860831969, +0.6612093864662645, +0.0000000000000000,
+0.6612093864662645, +0.6612093864662645, +0.0000000000000000,
+0.9324695142031521, +0.6612093864662645, +0.0000000000000000,
-0.9324695142031521, +0.9324695142031521, +0.0000000000000000,
-0.6612093864662645, +0.9324695142031521, +0.0000000000000000,
-0.2386191860831969, +0.9324695142031521, +0.0000000000000000,
+0.2386191860831969, +0.9324695142031521, +0.0000000000000000,
+0.6612093864662645, +0.9324695142031521, +0.0000000000000000,
+0.9324695142031521, +0.9324695142031521, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.0293520816889804,
+0.0618072933723834,
+0.0801651173178066,
+0.0801651173178066,
+0.0618072933723834,
+0.0293520816889804,
+0.0618072933723834,
+0.1301489125881675,
+0.1688053670875879,
+0.1688053670875879,
+0.1301489125881675,
+0.0618072933723834,
+0.0801651173178066,
+0.1688053670875879,
+0.2189434501672966,
+0.2189434501672966,
+0.1688053670875879,
+0.0801651173178066,
+0.0801651173178066,
+0.1688053670875879,
+0.2189434501672966,
+0.2189434501672966,
+0.1688053670875879,
+0.0801651173178066,
+0.0618072933723834,
+0.1301489125881675,
+0.1688053670875879,
+0.1688053670875879,
+0.1301489125881675,
+0.0618072933723834,
+0.0293520816889804,
+0.0618072933723834,
+0.0801651173178066,
+0.0801651173178066,
+0.0618072933723834,
+0.0293520816889804,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Quadrangle) && t_quadrature.hasOrd(7))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 49;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.9491079123427585, -0.9491079123427585, +0.0000000000000000,
-0.7415311855993945, -0.9491079123427585, +0.0000000000000000,
-0.4058451513773972, -0.9491079123427585, +0.0000000000000000,
+0.0000000000000000, -0.9491079123427585, +0.0000000000000000,
+0.4058451513773972, -0.9491079123427585, +0.0000000000000000,
+0.7415311855993945, -0.9491079123427585, +0.0000000000000000,
+0.9491079123427585, -0.9491079123427585, +0.0000000000000000,
-0.9491079123427585, -0.7415311855993945, +0.0000000000000000,
-0.7415311855993945, -0.7415311855993945, +0.0000000000000000,
-0.4058451513773972, -0.7415311855993945, +0.0000000000000000,
+0.0000000000000000, -0.7415311855993945, +0.0000000000000000,
+0.4058451513773972, -0.7415311855993945, +0.0000000000000000,
+0.7415311855993945, -0.7415311855993945, +0.0000000000000000,
+0.9491079123427585, -0.7415311855993945, +0.0000000000000000,
-0.9491079123427585, -0.4058451513773972, +0.0000000000000000,
-0.7415311855993945, -0.4058451513773972, +0.0000000000000000,
-0.4058451513773972, -0.4058451513773972, +0.0000000000000000,
+0.0000000000000000, -0.4058451513773972, +0.0000000000000000,
+0.4058451513773972, -0.4058451513773972, +0.0000000000000000,
+0.7415311855993945, -0.4058451513773972, +0.0000000000000000,
+0.9491079123427585, -0.4058451513773972, +0.0000000000000000,
-0.9491079123427585, +0.0000000000000000, +0.0000000000000000,
-0.7415311855993945, +0.0000000000000000, +0.0000000000000000,
-0.4058451513773972, +0.0000000000000000, +0.0000000000000000,
+0.0000000000000000, +0.0000000000000000, +0.0000000000000000,
+0.4058451513773972, +0.0000000000000000, +0.0000000000000000,
+0.7415311855993945, +0.0000000000000000, +0.0000000000000000,
+0.9491079123427585, +0.0000000000000000, +0.0000000000000000,
-0.9491079123427585, +0.4058451513773972, +0.0000000000000000,
-0.7415311855993945, +0.4058451513773972, +0.0000000000000000,
-0.4058451513773972, +0.4058451513773972, +0.0000000000000000,
+0.0000000000000000, +0.4058451513773972, +0.0000000000000000,
+0.4058451513773972, +0.4058451513773972, +0.0000000000000000,
+0.7415311855993945, +0.4058451513773972, +0.0000000000000000,
+0.9491079123427585, +0.4058451513773972, +0.0000000000000000,
-0.9491079123427585, +0.7415311855993945, +0.0000000000000000,
-0.7415311855993945, +0.7415311855993945, +0.0000000000000000,
-0.4058451513773972, +0.7415311855993945, +0.0000000000000000,
+0.0000000000000000, +0.7415311855993945, +0.0000000000000000,
+0.4058451513773972, +0.7415311855993945, +0.0000000000000000,
+0.7415311855993945, +0.7415311855993945, +0.0000000000000000,
+0.9491079123427585, +0.7415311855993945, +0.0000000000000000,
-0.9491079123427585, +0.9491079123427585, +0.0000000000000000,
-0.7415311855993945, +0.9491079123427585, +0.0000000000000000,
-0.4058451513773972, +0.9491079123427585, +0.0000000000000000,
+0.0000000000000000, +0.9491079123427585, +0.0000000000000000,
+0.4058451513773972, +0.9491079123427585, +0.0000000000000000,
+0.7415311855993945, +0.9491079123427585, +0.0000000000000000,
+0.9491079123427585, +0.9491079123427585, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.0167663564637535,
+0.0362176431542396,
+0.0494412511719133,
+0.0541194307579278,
+0.0494412511719133,
+0.0362176431542396,
+0.0167663564637535,
+0.0362176431542396,
+0.0782351060281695,
+0.1067999237589045,
+0.1169054370959262,
+0.1067999237589045,
+0.0782351060281695,
+0.0362176431542396,
+0.0494412511719133,
+0.1067999237589045,
+0.1457941874687416,
+0.1595893762111190,
+0.1457941874687416,
+0.1067999237589045,
+0.0494412511719133,
+0.0541194307579278,
+0.1169054370959262,
+0.1595893762111190,
+0.1746898792169928,
+0.1595893762111190,
+0.1169054370959262,
+0.0541194307579278,
+0.0494412511719133,
+0.1067999237589045,
+0.1457941874687416,
+0.1595893762111190,
+0.1457941874687416,
+0.1067999237589045,
+0.0494412511719133,
+0.0362176431542396,
+0.0782351060281695,
+0.1067999237589045,
+0.1169054370959262,
+0.1067999237589045,
+0.0782351060281695,
+0.0362176431542396,
+0.0167663564637535,
+0.0362176431542396,
+0.0494412511719133,
+0.0541194307579278,
+0.0494412511719133,
+0.0362176431542396,
+0.0167663564637535,
        };

    };

    template<Element t_element, Quadrature t_quadrature>
    requires(t_element.hasShape(Element::Shape::Quadrangle) && t_quadrature.hasOrd(8))
    struct ElementQuadratureRuleTraits<t_element, t_quadrature>
    {

        static constexpr
        Integer
        getSize()
        {
            return 64;
        }

        static constexpr
        std::array<std::array<Real, 3>, getSize()> const &
        getReferencePoints()
        {
            return reference_points_;
        }

        static constexpr
        std::array<Real, getSize()> const &
        getReferenceWeights()
        {
            return reference_weights_;
        }

        std::array<std::array<Real, 3>, getSize()> static constexpr reference_points_ = {
-0.9602898564975362, -0.9602898564975362, +0.0000000000000000,
-0.7966664774136267, -0.9602898564975362, +0.0000000000000000,
-0.5255324099163290, -0.9602898564975362, +0.0000000000000000,
-0.1834346424956498, -0.9602898564975362, +0.0000000000000000,
+0.1834346424956498, -0.9602898564975362, +0.0000000000000000,
+0.5255324099163290, -0.9602898564975362, +0.0000000000000000,
+0.7966664774136267, -0.9602898564975362, +0.0000000000000000,
+0.9602898564975362, -0.9602898564975362, +0.0000000000000000,
-0.9602898564975362, -0.7966664774136267, +0.0000000000000000,
-0.7966664774136267, -0.7966664774136267, +0.0000000000000000,
-0.5255324099163290, -0.7966664774136267, +0.0000000000000000,
-0.1834346424956498, -0.7966664774136267, +0.0000000000000000,
+0.1834346424956498, -0.7966664774136267, +0.0000000000000000,
+0.5255324099163290, -0.7966664774136267, +0.0000000000000000,
+0.7966664774136267, -0.7966664774136267, +0.0000000000000000,
+0.9602898564975362, -0.7966664774136267, +0.0000000000000000,
-0.9602898564975362, -0.5255324099163290, +0.0000000000000000,
-0.7966664774136267, -0.5255324099163290, +0.0000000000000000,
-0.5255324099163290, -0.5255324099163290, +0.0000000000000000,
-0.1834346424956498, -0.5255324099163290, +0.0000000000000000,
+0.1834346424956498, -0.5255324099163290, +0.0000000000000000,
+0.5255324099163290, -0.5255324099163290, +0.0000000000000000,
+0.7966664774136267, -0.5255324099163290, +0.0000000000000000,
+0.9602898564975362, -0.5255324099163290, +0.0000000000000000,
-0.9602898564975362, -0.1834346424956498, +0.0000000000000000,
-0.7966664774136267, -0.1834346424956498, +0.0000000000000000,
-0.5255324099163290, -0.1834346424956498, +0.0000000000000000,
-0.1834346424956498, -0.1834346424956498, +0.0000000000000000,
+0.1834346424956498, -0.1834346424956498, +0.0000000000000000,
+0.5255324099163290, -0.1834346424956498, +0.0000000000000000,
+0.7966664774136267, -0.1834346424956498, +0.0000000000000000,
+0.9602898564975362, -0.1834346424956498, +0.0000000000000000,
-0.9602898564975362, +0.1834346424956498, +0.0000000000000000,
-0.7966664774136267, +0.1834346424956498, +0.0000000000000000,
-0.5255324099163290, +0.1834346424956498, +0.0000000000000000,
-0.1834346424956498, +0.1834346424956498, +0.0000000000000000,
+0.1834346424956498, +0.1834346424956498, +0.0000000000000000,
+0.5255324099163290, +0.1834346424956498, +0.0000000000000000,
+0.7966664774136267, +0.1834346424956498, +0.0000000000000000,
+0.9602898564975362, +0.1834346424956498, +0.0000000000000000,
-0.9602898564975362, +0.5255324099163290, +0.0000000000000000,
-0.7966664774136267, +0.5255324099163290, +0.0000000000000000,
-0.5255324099163290, +0.5255324099163290, +0.0000000000000000,
-0.1834346424956498, +0.5255324099163290, +0.0000000000000000,
+0.1834346424956498, +0.5255324099163290, +0.0000000000000000,
+0.5255324099163290, +0.5255324099163290, +0.0000000000000000,
+0.7966664774136267, +0.5255324099163290, +0.0000000000000000,
+0.9602898564975362, +0.5255324099163290, +0.0000000000000000,
-0.9602898564975362, +0.7966664774136267, +0.0000000000000000,
-0.7966664774136267, +0.7966664774136267, +0.0000000000000000,
-0.5255324099163290, +0.7966664774136267, +0.0000000000000000,
-0.1834346424956498, +0.7966664774136267, +0.0000000000000000,
+0.1834346424956498, +0.7966664774136267, +0.0000000000000000,
+0.5255324099163290, +0.7966664774136267, +0.0000000000000000,
+0.7966664774136267, +0.7966664774136267, +0.0000000000000000,
+0.9602898564975362, +0.7966664774136267, +0.0000000000000000,
-0.9602898564975362, +0.9602898564975362, +0.0000000000000000,
-0.7966664774136267, +0.9602898564975362, +0.0000000000000000,
-0.5255324099163290, +0.9602898564975362, +0.0000000000000000,
-0.1834346424956498, +0.9602898564975362, +0.0000000000000000,
+0.1834346424956498, +0.9602898564975362, +0.0000000000000000,
+0.5255324099163290, +0.9602898564975362, +0.0000000000000000,
+0.7966664774136267, +0.9602898564975362, +0.0000000000000000,
+0.9602898564975362, +0.9602898564975362, +0.0000000000000000,
        };
        
        std::array<Real, getSize()> static constexpr reference_weights_ = {
+0.0102472165594920,
+0.0225113066164548,
+0.0317560645867820,
+0.0367139485276475,
+0.0367139485276475,
+0.0317560645867820,
+0.0225113066164548,
+0.0102472165594920,
+0.0225113066164548,
+0.0494533244845528,
+0.0697624084252229,
+0.0806539949271436,
+0.0806539949271436,
+0.0697624084252229,
+0.0494533244845528,
+0.0225113066164548,
+0.0317560645867820,
+0.0697624084252229,
+0.0984118596679543,
+0.1137763131979283,
+0.1137763131979283,
+0.0984118596679543,
+0.0697624084252229,
+0.0317560645867820,
+0.0367139485276475,
+0.0806539949271436,
+0.1137763131979283,
+0.1315395267256426,
+0.1315395267256426,
+0.1137763131979283,
+0.0806539949271436,
+0.0367139485276475,
+0.0367139485276475,
+0.0806539949271436,
+0.1137763131979283,
+0.1315395267256426,
+0.1315395267256426,
+0.1137763131979283,
+0.0806539949271436,
+0.0367139485276475,
+0.0317560645867820,
+0.0697624084252229,
+0.0984118596679543,
+0.1137763131979283,
+0.1137763131979283,
+0.0984118596679543,
+0.0697624084252229,
+0.0317560645867820,
+0.0225113066164548,
+0.0494533244845528,
+0.0697624084252229,
+0.0806539949271436,
+0.0806539949271436,
+0.0697624084252229,
+0.0494533244845528,
+0.0225113066164548,
+0.0102472165594920,
+0.0225113066164548,
+0.0317560645867820,
+0.0367139485276475,
+0.0367139485276475,
+0.0317560645867820,
+0.0225113066164548,
+0.0102472165594920,
        };

    };

} // namespace lolita

#endif /* CC79BDC5_49DB_4A81_8A93_18ABD6551AF1 */
