#ifndef SCENE_2D
#define SCENE_2D

#include "Scene2d.hpp"
#include <iostream>

using namespace CPGF::Objects2d;
using namespace CPGF;

Scene2d::Scene2d(const std::vector<Object2d*>& objects,
            const std::vector<Text*>& texts):
    objects(objects), texts(texts)
{

}

Scene2d& Scene2d::add(Object2d& object)
{
    objects.push_back(&object);
    return *this;
}

Scene2d& Scene2d::add(Text& text)
{
    texts.push_back(&text);
    return *this;
}

Scene2d operator+(const Scene2d& S1, const Scene2d& S2)
{
    Scene2d S = S1;
    S += S2;
    return S;
}

Scene2d& Scene2d::operator+=(const Scene2d& S2)
{
    for (unsigned int i = 0; i < S2.objects.size(); i++)
    {
        objects.push_back(S2.objects[i]);
    }
    for (unsigned int i = 0; i < S2.texts.size(); i++)
    {
        texts.push_back(S2.texts[i]);
    }
    return *this;
}

Scene2d& Scene2d::operator+=(Object2d& Obj)
{
    add(Obj);
    return *this;
}

Scene2d& Scene2d::operator+=(Text& text)
{
    add(text);
    return *this;
}

std::string Scene2d::render_to_string() const
{
    std::string text = "\\documentclass[margin=4mm]{standalone}\n\\usepackage{pgf}\n\\usepackage{amssymb}\n\\usepackage{amsmath}\n\\begin{document}\n\\begin{pgfpicture}\n";

    for (unsigned int i = 0; i < objects.size(); i++)
    {
        text += objects[i]->render_to_string();
    }

    for (unsigned int i = 0; i < texts.size(); i++)
    {
        text += texts[i]->render_to_string();
    }

    text += "\\end{pgfpicture}\n\\end{document}\n";
    return text;
}

void Scene2d::render(const std::string& filename, const unsigned int density) const
{
    unsigned int dot_pos = filename.find_last_of(".");
    std::string extension = filename.substr(dot_pos + 1);
    unsigned int slash_pos = filename.find_last_of("/");
    std::string name;
    std::string directory;
    std::string pre_command;
    if (slash_pos > filename.length())
    {
        name = filename.substr(0, dot_pos);
        directory = ".";
        pre_command = "";
    }
    else
    {
        name = filename.substr(slash_pos+1, dot_pos - slash_pos - 1);
        directory = filename.substr(0, slash_pos);
        pre_command = "cd " + directory + "; ";
    }

    // We write the text in a file.
    FILE* F = fopen((directory + "/" + name + ".tex").c_str(), "w");
    fprintf(F, "%s", render_to_string().c_str());
    fclose(F);

    // We compile the latex file.
    system(( pre_command
        + "lualatex -output-format='pdf' " + name + ".tex >>/dev/null").c_str());
    // We eliminate unwanted files.
    system(( pre_command
        + "rm -r " + name + ".aux " + name + ".log >>/dev/null").c_str());

    if (extension != "pdf")
    {
        system(( pre_command
            + "convert -density " + std::to_string(density) + " " + name + ".pdf " + name + "." + extension
            + ">>/dev/null").c_str());
        system(( pre_command
            + "rm -r " + name + ".pdf >>/dev/null").c_str());
    }

    system(( pre_command
        + "rm -r " + name + ".tex >>/dev/null").c_str());

}

#endif // SCENE_2D
