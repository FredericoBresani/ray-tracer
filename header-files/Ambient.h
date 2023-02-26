#ifndef __AMBIENT__
#define __AMBIENT__

#include "Vectors.h"

class Ambient 
{
    public:
        Vec3D color;
        Ambient(const Vec3D &c): color(c) {}
        ~Ambient() {}
};

#endif