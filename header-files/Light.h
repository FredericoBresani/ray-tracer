#ifndef __LIGHT__
#define __LIGHT__

#include "Vectors.h"
#include "Points.h"

class Light
{
    public:
        Point3D lightPos;
        Vec3D lightColor;
        Light(const Point3D& pos, const Vec3D& color): lightPos(pos), lightColor(color) {}
        ~Light() {}
};

#endif