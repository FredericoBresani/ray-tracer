#ifndef __AMBIENT__
#define __AMBIENT__

#include "RGBColor.h"

class Ambient 
{
    public:
        RGBColor color;
        Ambient(const RGBColor &c): color(c) {}
        ~Ambient() {}
};

#endif