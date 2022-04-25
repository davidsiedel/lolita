//
// Created by dsiedel on 17/04/22.
//

#ifndef LOLITA_LOLITA_FILE_HXX
#define LOLITA_LOLITA_FILE_HXX

#include <fstream>

#include "lolita/lolita.hxx"
#include "lolita/lolita_lolita.hxx"
#include "lolita/lolita_containers.hxx"

namespace lolita::file
{

    static void inline
    removeCharacter(
            Strg &
            line,
            Char &&
            c
    )
    {
        line.erase(std::remove(line.begin(), line.end(), c), line.end());
    }

    struct File : public Array<Strg>
    {

        File() = default;

        explicit
        File(
                Strg const &
                file_path_arg
        )
        :
        Array<Strg>(readLines(file_path_arg))
        {}

        Bool
        operator==(
                File const &
                other
        )
        const = default;

        Bool
        operator!=(
                File const &
                other
        )
        const = default;

    private:

        static inline
        Array<Strg>
        readLines(
                Strg const &
                f
        )
        {
            FileStream file(f);
            if (!file) {
                throw std::runtime_error("Could not open file");
            }
            Array<Strg> lines;
            for (Strg line; std::getline(file, line); ) {
                lines.data.push_back(line);
            }
            return lines;
        }

    };

    using File = file::File;

}

#endif //LOLITA_LOLITA_FILE_HXX
