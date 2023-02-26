#ifndef __LIGHT__
#define __LIGHT__

#include "Vectors.h"
#include "RGBColor.h"

class Light
{
    public:
        Point3D lightPos;
        RGBColor lightColor;
        Light(const Point3D& pos, const RGBColor& color): lightPos(pos), lightColor(color) {}
        ~Light() {}
};

#endif