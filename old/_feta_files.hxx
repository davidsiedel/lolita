//
// Created by dsiedel on 20/11/2021.
//

#ifndef FETA__FETA_FILES_HXX
#define FETA__FETA_FILES_HXX

#include <fstream>
#include "new/_feta_array.hxx"
#include "new/_feta.hxx"
#include "new/_feta_raise.hxx"

namespace feta::file
{

    enum struct FileType
    {

        Input,
        Output

    };

    void inline
    removeCharacter(
            Strg &
            line,
            Char &&
            c
    )
    {
        line.erase(std::remove(line.begin(), line.end(), c), line.end());
    }

    template<FileType>
    struct File;

    template<>
    struct File<FileType::Input>
    {

        File() = default;

        explicit File(
                const Strg &
                f
        )
        :
        content(readLines(f))
        {}

        Bool
        operator==(
                File const &
                other
        )
        const
        {
            auto eq_0 = file_path == other.file_path;
            return eq_0;
        }

        Bool
        operator!=(
                File const &
                other
        )
        const
        {
            return !(other == * this);
        }

        Indx
        getLineIndex(
                Strg const &
                line
        )
        const
        {
            Indx line_index = content.index(line);
            if (line_index == content.getSize()) {
                aassert(true, "line was not found");
            }
            return line_index;
        }

        Strg const &
        getLine(
                Indx
                i
        )
        const
        {
            return content(i);
        }

        Strg &
        getLine(
                Indx
                i
        )
        {
            return content(i);
        }

        Strg const &
        getFilePath()
        const
        {
            return file_path;
        }

    private:

        static
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

        Strg file_path;

        Array<Strg> content;

    };

}

#endif //FETA__FETA_FILES_HXX
