#ifndef __AMBIENT__
#define __AMBIENT__

#include "RGBColor.h"

class Ambient 
{
    public:
        RGBColor color;
        float ir;
        Ambient(const RGBColor &c, const float &r): color(c), ir(r) {}
        ~Ambient() {}
};

#endif