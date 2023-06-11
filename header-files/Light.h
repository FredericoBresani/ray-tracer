#ifndef __LIGHT__
#define __LIGHT__

#include "Vectors.h"
#include "RGBColor.h"

class Light
{
    public:
        //shadows
        Point3D lightPos;
        RGBColor lightColor;


        virtual Vec3D getDirection() = 0;
        virtual RGBColor incidentRadiance() = 0;
};

#endif