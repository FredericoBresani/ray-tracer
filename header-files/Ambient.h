#ifndef __AMBIENT__
#define __AMBIENT__

#include "RGBColor.h"

class Ambient 
{
    public:
        RGBColor color;
        float ir;
        int depth = 1;
        Ambient(const RGBColor &c, const float &r, const int &d): color(c), ir(r), depth(d) {}
        ~Ambient() {}
};

#endif