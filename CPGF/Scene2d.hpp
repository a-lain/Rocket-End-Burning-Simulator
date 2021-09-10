#ifndef SCENE2D_HPP
#define SCENE2D_HPP

#include <vector>
#include <string>

#include "Objects2d/Object2d.hpp"
#include "Text/Text.hpp"

namespace CPGF
{
    class Scene2d
    {
        public:
        std::vector<Objects2d::Object2d*> objects;
        std::vector<Text*> texts;

        Scene2d(const std::vector<Objects2d::Object2d*>& objects = std::vector<Objects2d::Object2d*>(),
            const std::vector<Text*>& texts = std::vector<Text*>());

        Scene2d& add(Objects2d::Object2d& object);
        Scene2d& add(Text& texts);

        Scene2d friend operator+(const Scene2d& S1, const Scene2d& S2);
        Scene2d& operator+=(const Scene2d& S2);
        Scene2d& operator+=(Objects2d::Object2d& Obj);
        Scene2d& operator+=(Text& text);

        std::string render_to_string() const;
        void render(const std::string& filename, const unsigned int density = 100) const;
    };
}

#endif // SCENE2D_HPP