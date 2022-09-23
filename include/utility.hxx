#ifndef BD1D6E93_67D6_418F_A38E_F9C11C70DB4B
#define BD1D6E93_67D6_418F_A38E_F9C11C70DB4B

#include <iostream>
#include <fstream>
#include <execution>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <map>
#include <ostream>
#include <iomanip>
#include <filesystem>
#include <span>

#include "config.hxx"

namespace lolita::utility
{

    template<typename t_T>
    static
    std::vector<t_T>
    linspace(
        t_T start,
        t_T end,
        Integer steps
    )
    {
        auto linspacee = std::vector<t_T>(steps);
        auto step = (end - start) / (steps - 1);
        auto count = 0;
        for (auto & val : linspacee)
        {
            val = start + count * step;
            count ++;
        }
        return linspacee;
    }

    struct TimeSteps
    {

        TimeSteps()
        :
        step_(0),
        iteration_(0),
        values_(1, 0)
        {}

        inline
        void
        addRange(
            Real end,
            Integer steps
        )
        {
            auto start = values_.back();
            if (end > start)
            {
                auto step = (end - start) / (steps);
                for (auto i = 1; i < steps + 1; i++)
                {
                    values_.push_back(start + i * step);
                }
            }
            else
            {
                throw std::runtime_error("NO !");
            }
        }
        
        inline
        Real const &
        getTime()
        const
        {
            return values_[step_];
        }
        
        inline
        Integer
        getStep()
        const
        {
            return step_;
        }
        
        inline
        Integer
        getIteration()
        const
        {
            return iteration_;
        }
        
        inline
        Boolean
        isDone()
        const
        {
            return step_ == values_.size() - 1;
        }
        
        inline
        Boolean
        hasIteration(
            Integer iteration
        )
        const
        {
            return iteration_ == iteration;
        }

        inline
        void
        split()
        {
            if (step_ > 0)
            {
                auto value = values_[step_ - 1] + (1.0 / 2.0) * (values_[step_] - values_[step_ - 1]);
                values_.insert(values_.begin() + step_, value);
                iteration_ = 0;
            }
            else
            {
                throw std::runtime_error("NO !");
            }
        }

        inline
        void
        update()
        {
            step_ ++;
            iteration_ = 0;
        }
        
        inline
        void
        updateIteration()
        {
            iteration_ ++;
        }
        
        inline
        friend
        std::ostream &
        operator<<(
            std::ostream & os,
            TimeSteps const & time_steps
        )
        {
            os << "time_steps : ";
            for (auto time_step : time_steps.values_)
            {
                os << time_step << " / ";
            }
            return os;
        }

    private:

        Integer step_;

        Integer iteration_;

        std::vector<Real> values_;

    };

    template<typename t_T>
    struct Holderr
    {

    private:

        struct Item
        {

            template<typename... t_U>
            static
            Item
            make(
                std::basic_string_view<Character> label,
                t_U const &... args
            )
            {
                return Item(label, t_T(args...));
            }

            template<typename... t_U>
            static
            Item
            make(
                std::basic_string_view<Character> label,
                t_U &&... args
            )
            {
                return Item(label, t_T(std::forward<t_U>(args)...));
            }

            Item(
                std::basic_string_view<Character> label,
                t_T const & value
            )
            :
            label_(label),
            value_(value)
            {}

            Item(
                std::basic_string_view<Character> label,
                t_T && value
            )
            :
            label_(label),
            value_(std::forward<t_T>(value))
            {}

            t_T const &
            getValue()
            const
            {
                return value_;
            }

            t_T &
            getValue()
            {
                return value_;
            }

            std::basic_string_view<Character>
            getLabel()
            const
            {
                return label_;
            }

            t_T value_;

            std::basic_string_view<Character> label_;

        };

        Integer
        getItemIndex(
            std::basic_string_view<Character> label
        )
        const
        {
            return std::distance(items_.begin(), std::find_if(items_.begin(), items_.end(), [&] (Item const & item) { return item.getLabel() == label; }));
        }

    public:

        Holderr()
        {}

        t_T const &
        getValue(
            std::basic_string_view<Character> label
        )
        const
        {
            auto item_index = getItemIndex(label);
            if (item_index != items_.size())
            {
                return items_[item_index].getValue();
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }

        t_T &
        getValue(
            std::basic_string_view<Character> label
        )
        {
            auto item_index = getItemIndex(label);
            if (item_index != items_.size())
            {
                return items_[item_index].getValue();
            }
            else
            {
                throw std::runtime_error("NO");
            }
        }

        // void
        // setItem(
        //     std::basic_string_view<Character> label,
        //     t_T const & value
        // )
        // {
        //     auto item_index = getItemIndex(label);
        //     if (item_index != items_.size())
        //     {
        //         items_[item_index] = Item(label, value);
        //     }
        //     else
        //     {
        //         items_.push_back(Item(label, value));
        //     }
        // }

        // void
        // setItem(
        //     std::basic_string_view<Character> label,
        //     t_T && value
        // )
        // {
        //     auto item_index = getItemIndex(label);
        //     if (item_index != items_.size())
        //     {
        //         items_[item_index] = Item(label, std::forward<t_T>(value));
        //     }
        //     else
        //     {
        //         items_.push_back(Item(label, std::forward<t_T>(value)));
        //     }
        // }

        template<typename... t_U>
        void
        setItem(
            std::basic_string_view<Character> label,
            t_U const &... args
        )
        {
            auto item_index = getItemIndex(label);
            if (item_index != items_.size())
            {
                items_[item_index] = Item::template make<t_U...>(label, args...);
            }
            else
            {
                items_.push_back(Item::template make<t_U...>(label, args...));
            }
        }

        template<typename... t_U>
        void
        setItem(
            std::basic_string_view<Character> label,
            t_U &&... args
        )
        {
            auto item_index = getItemIndex(label);
            if (item_index != items_.size())
            {
                items_[item_index] = Item::template make<t_U...>(label, std::forward<t_U>(args)...);
            }
            else
            {
                items_.push_back(Item::template make<t_U...>(label, std::forward<t_U>(args)...));
            }
        }

    private:

        std::vector<Item> items_;

    };

    template<auto t_value>
    struct Holder
    {

        auto static constexpr value_ = t_value;

    };

    template<typename t_T, typename t_U>
    static constexpr
    Boolean
    areEqual(
        t_T && x,
        t_U && y
    )
    {
        if constexpr (std::convertible_to<t_T, t_U>)
        {
            if (std::forward<t_T>(x) == std::forward<t_U>(y))
            {
                return true;
            }
            else
            {
                return false;
            }
        }
        else
        {
            return false;
        }
    }
    
    template<typename t_T, typename... t_U>
    static constexpr
    Boolean
    isAnyOf(
        t_T && x,
        t_U &&... y
    )
    {
        return ((areEqual(std::forward<t_T>(x), std::forward<t_U>(y)) || ...));
    }

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, typename t_U>
    struct tuple_concatenation_traits;

    template<typename... t_T, typename... t_U>
    struct tuple_concatenation_traits<std::tuple<t_T...>, std::tuple<t_U...>>
    {

        using type = std::tuple<t_T..., t_U...>;

    };

    template<typename t_T, typename t_U>
    using tuple_cat_t = typename tuple_concatenation_traits<t_T, t_U>::type;

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, Integer t_a, Integer t_b>
    struct TupleSliceTraits;

    template<typename... t_T, Integer t_a, Integer t_b>
    struct TupleSliceTraits<std::tuple<t_T...>, t_a, t_b>
    {

        using type = tuple_cat_t<std::tuple<std::tuple_element_t<t_a, std::tuple<t_T...>>>, typename TupleSliceTraits<std::tuple<t_T...>, t_a + 1, t_b>::type>;

    };
    
    template<typename... t_T, Integer t_a, Integer t_b>
    requires(t_a == t_b - 1)
    struct TupleSliceTraits<std::tuple<t_T...>, t_a, t_b>
    {

        using type = std::tuple<std::tuple_element_t<t_b - 1, std::tuple<t_T...>>>;

    };
    
    template<typename... t_T, Integer t_a, Integer t_b>
    requires(t_a == t_b)
    struct TupleSliceTraits<std::tuple<t_T...>, t_a, t_b>
    {

        using type = std::tuple<>;

    };

    template<typename t_T, Integer t_a, Integer t_b>
    using tuple_slice_t = typename TupleSliceTraits<t_T, t_a, t_b>::type;

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, typename t_U>
    struct TupleHasTraits;

    template<typename... t_T, typename t_U>
    struct TupleHasTraits<std::tuple<t_T...>, t_U>
    {

        Boolean static constexpr value = (std::is_same_v<t_U, t_T> || ...);

    };

    template<typename t_T, typename t_U>
    static constexpr
    Boolean tuple_has_v = TupleHasTraits<t_T, t_U>::value;

    // ------------------------------------------------------------------------------------------------------

    template<typename t_T, Integer t_a>
    struct TupleUniqueTraits;

    template<typename... t_T, Integer t_a>
    requires(t_a < sizeof...(t_T) - 1 && !tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = tuple_cat_t<std::tuple<std::tuple_element_t<t_a, std::tuple<t_T...>>>, typename TupleUniqueTraits<std::tuple<t_T...>, t_a + 1>::type>;

    };

    template<typename... t_T, Integer t_a>
    requires(t_a < sizeof...(t_T) - 1 && tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = tuple_cat_t<std::tuple<>, typename TupleUniqueTraits<std::tuple<t_T...>, t_a + 1>::type>;

    };

    template<typename... t_T, Integer t_a>
    requires(t_a == sizeof...(t_T) - 1 && !tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = std::tuple<std::tuple_element_t<t_a, std::tuple<t_T...>>>;

    };

    template<typename... t_T, Integer t_a>
    requires(t_a == sizeof...(t_T) - 1 && tuple_has_v<tuple_slice_t<std::tuple<t_T...>, 0, t_a>, std::tuple_element_t<t_a, std::tuple<t_T...>>>)
    struct TupleUniqueTraits<std::tuple<t_T...>, t_a>
    {

        using type = std::tuple<>;

    };

    template<typename t_T>
    using tuple_unique_t = typename TupleUniqueTraits<t_T, 0>::type;

    // ------------------------------------------------------------------------------------------------------

    template<typename>
    void
    TD()
    {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
    }

    // ------------------------------------------------------------------------------------------------------

    template<Integer t_i, typename t_T>
    struct AggregateTail
    {

        constexpr
        AggregateTail()
        :
        value_(t_T())
        {}

        constexpr explicit
        AggregateTail(
            t_T const & value
        )
        :
        value_(value)
        {}

        constexpr explicit
        AggregateTail(
            t_T && value
        )
        :
        value_(std::move(value))
        {}

        constexpr
        Boolean
        operator==(
            AggregateTail const &
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            AggregateTail const &
        )
        const = default;

        t_T value_;

    };

    template<Integer t_i, typename... t_T>
    struct AggregateHead;

    template<Integer t_i>
    struct AggregateHead<t_i>
    {

        constexpr
        AggregateHead()
        {}

        constexpr
        Boolean
        operator==(
            AggregateHead const &
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            AggregateHead const &
        )
        const = default;

    };

    template<Integer t_i, typename t_T, typename... t_U>
    struct AggregateHead<t_i, t_T, t_U...> : AggregateTail<t_i, t_T>, AggregateHead<t_i + 1, t_U...>
    {

        constexpr
        AggregateHead(
            t_T const & value,
            t_U const &... values
        )
        :
        AggregateTail<t_i, t_T>(value),
        AggregateHead<t_i + 1, t_U...>(values...)
        {}

        constexpr
        AggregateHead(
            t_T && value,
            t_U &&... values
        )
        :
        AggregateTail<t_i, t_T>(std::move(value)),
        AggregateHead<t_i + 1, t_U...>(std::move(values)...)
        {}

        constexpr
        Boolean
        operator==(
            AggregateHead const &
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            AggregateHead const &
        )
        const = default;

    };

    template<Integer t_i, typename t_T, typename... t_U>
    static constexpr
    t_T &
    get(
        AggregateHead<t_i, t_T, t_U...> & aggregate
    )
    {
        return aggregate.AggregateTail<t_i, t_T>::value_;
    }

    template<Integer t_i, typename t_T, typename... t_U>
    static constexpr
    t_T const &
    get(
        AggregateHead<t_i, t_T, t_U...> const & aggregate
    )
    {
        return aggregate.AggregateTail<t_i, t_T>::value_;
    }

    template<typename... T>
    using Aggregate = AggregateHead<0, T...>;

    // ------------------------------------------------------------------------------------------------------

    template<typename T>
    struct AggregateSizeTraits;

    template<typename... T>
    struct AggregateSizeTraits<Aggregate<T...>>
    {

        Integer static constexpr value = sizeof...(T);

    };

    template<typename T>
    static constexpr
    Integer aggregate_size_v = AggregateSizeTraits<T>::value;

    template<Integer I, typename T>
    struct AggregateElementTraits;

    template<Integer I, typename... T>
    struct AggregateElementTraits<I, Aggregate<T...>>
    {

        using type = std::tuple_element_t<I, std::tuple<T...>>;

    };

    template<Integer I, typename... T>
    using aggregate_element_t = typename AggregateElementTraits<I, T...>::type;

    template<template<auto> typename t_T, auto t_aggregate, Integer t_i>
    struct aggregate_expansion_traits
    {

        using type = tuple_cat_t<std::tuple<t_T<get<t_i>(t_aggregate)>>, typename aggregate_expansion_traits<t_T, t_aggregate, t_i + 1>::type>;

    };
    
    template<template<auto> typename t_T, auto t_aggregate, Integer t_i>
    requires(t_i == aggregate_size_v<std::decay_t<decltype(t_aggregate)>> - 1)
    struct aggregate_expansion_traits<t_T, t_aggregate, t_i>
    {

        using type = std::tuple<t_T<get<t_i>(t_aggregate)>>;

    };

    template<template<auto> typename t_T, auto t_aggregate>
    using aggregate_expansion_t = typename aggregate_expansion_traits<t_T, t_aggregate, 0>::type;

    // ------------------------------------------------------------------------------------------------------

    struct Label
    {

    private:

        Integer static constexpr size_ = 50;

        static constexpr
        std::array<Character, size_>
        tag(
            std::basic_string_view<Character> && str
        )
        {
            auto tag = std::array<Character, size_>();
            for (auto i = 0; i < std::forward<std::basic_string_view<Character>>(str).size(); i++)
            {
                tag[i] = std::forward<std::basic_string_view<Character>>(str)[i];
            }
            return tag;
        }

    public:

        constexpr
        Label()
        :
        tag_()
        {}

        explicit constexpr
        Label(
            std::basic_string_view<Character> && str
        )
        :
        tag_(tag(std::forward<std::basic_string_view<Character>>(str)))
        {}

        // constexpr
        // Label(
        //     Character const * str
        // )
        // :
        // tag_(tag(std::basic_string_view<Character>(str)))
        // {}

        constexpr
        Boolean
        operator==(
            Label const & other
        )
        const = default;

        constexpr
        Boolean
        operator!=(
            Label const & other
        )
        const = default;

        // constexpr
        // Boolean
        // operator==(
        //     auto const & str
        // )
        // const
        // {
        //     return str == this->view();
        // }

        // constexpr
        // Boolean
        // operator!=(
        //     auto const & str
        // )
        // const
        // {
        //     return !(* this == str);
        // }

        constexpr
        Boolean
        operator==(
            auto && str
        )
        const
        {
            return this->view() == std::forward<decltype(str)>(str);
        }

        constexpr
        Boolean
        operator!=(
            auto && str
        )
        const
        {
            return !(* this == str);
        }

        constexpr
        std::basic_string_view<Character>
        view()
        const
        {
            return std::basic_string_view<Character>(tag_.data(), std::distance(tag_.begin(), std::find(tag_.begin(), tag_.end(), Character())));
        }

        friend inline
        std::ostream &
        operator<<(
            std::ostream & os,
            Label const & label
        )
        {
            os << label.view();
            return os;
        }

        std::array<Character, size_> tag_;

    };

    template<typename... t_T>
    static constexpr
    std::array<Label, sizeof...(t_T)>
    makeLabels(
        t_T &&... strs
    )
    {
        return std::array<Label, sizeof...(t_T)>{Label(std::forward<std::basic_string_view<Character>>(strs))...};
    }

    // ------------------------------------------------------------------------------------------------------

    // struct File2
    // {

    //     static inline
    //     File2
    //     make(
    //         std::basic_string<Character> const & file_path
    //     )
    //     {
    //         auto file_stream = std::basic_ofstream<Character>();
    //         file_stream.open(file_path);
    //     }

    //     static inline
    //     File2
    //     open(
    //         std::basic_string<Character> const & file_path
    //     )
    //     {
    //         auto file_stream = std::basic_ofstream<Character>();
    //         file_stream.open(file_path, std::ios_base::app);
    //     }

    //     std::basic_ofstream<Character>;

    // };

    struct File
    {

        // static inline
        // File
        // make(
        //     std::basic_string_view<Character> file_path
        // )
        // {
            
        // }

        // static inline
        // File
        // open(
        //     std::basic_string<Character> const & file_path
        // )
        // {
        //     auto outfile = std::basic_ifstream<Character>();
        //     outfile.open(file_path);
        //     auto lines = 
        // }
        
        // inline
        // void
        // addLine()
        // {

        // }

        File() = default;

        explicit
        File(
            std::basic_string_view<Character> file_path
        )
        :
        file_path_(file_path)
        {
            readLines();
        }

        Boolean
        operator==(
            File const &
            other
        )
        const = default;

        Boolean
        operator!=(
            File const &
            other
        )
        const = default;

    private:

        inline
        void
        readLines()
        {
            std::basic_ifstream<Character> file(file_path_);
            if (!file)
            {
                throw std::runtime_error("Could not open file");
            }
            for (std::basic_string<Character> line; std::getline(file, line); )
            {
                lines_.push_back(line);
            }
        }

    public:

        std::vector<std::basic_string<Character>> lines_;

        std::basic_string<Character> file_path_;

    };

} // namespace lolita::utility


#endif /* BD1D6E93_67D6_418F_A38E_F9C11C70DB4B */
