#ifndef DF64FB37_2673_4868_AC2E_AAD4C50FECD0
#define DF64FB37_2673_4868_AC2E_AAD4C50FECD0

#include "config.hxx"

namespace lolita
{
    
    struct Label
    {

    private:

        Integer static constexpr size_ = 100;

        static constexpr
        std::array<Character, size_>
        tag(
            auto &&... str
        )
        {
            auto tag = std::array<Character, size_>();
            auto offset = Integer(0);
            auto make_tag = [&] (
                std::basic_string_view<Character> && str_arg
            )
            constexpr
            {
                for (auto i = 0; i < std::forward<std::basic_string_view<Character>>(str_arg).size(); i++)
                {
                    tag[offset + i] = std::forward<std::basic_string_view<Character>>(str_arg)[i];
                }
                offset += std::forward<std::basic_string_view<Character>>(str_arg).size();
            };
            ((make_tag(std::forward<std::basic_string_view<Character>>(str)), ...));
            return tag;
        }

    public:

        constexpr
        Label(
            auto &&... str
        )
        :
        tag_(tag(std::forward<std::basic_string_view<Character>>(str)...))
        {}

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
            return !(* this == std::forward<decltype(str)>(str));
        }

        constexpr
        Boolean
        is(
            std::basic_string_view<Character> && str
        )
        const
        {
            return * this == std::forward<std::basic_string_view<Character>>(str);
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

} // namespace lolita


#endif /* DF64FB37_2673_4868_AC2E_AAD4C50FECD0 */
