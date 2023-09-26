#ifndef __LIGHT__
#define __LIGHT__

#include "Vectors.h"
#include "RGBColor.h"
#include "HitInfo.h"
#include "Points.h"

class Light
{
    public:
        virtual Vec3D getDirection(HitInfo &hit) = 0;
        virtual RGBColor incidentRadiance(HitInfo &hit) = 0;
        virtual Point3D getPos() = 0;
        virtual RGBColor getColor() = 0;
        virtual bool castShadows() = 0;
    protected:
        bool shadows;
};

#endif