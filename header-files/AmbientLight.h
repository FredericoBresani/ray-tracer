#ifndef __AMBIENTLIGHT__
#define __AMBIENTLIGHT__


#include "Light.h"
#include "Points.h"
#include "Vectors.h"
#include "HitInfo.h"

class AmbientLight: public Light {

    public:
        Vec3D getDirection(HitInfo &hit);
        RGBColor incidentRadiance(HitInfo &hit);
        Point3D getPos();
        RGBColor getColor();
    private:
        RGBColor lightColor;
        double reflectiveK;
};

Vec3D AmbientLight::getDirection(HitInfo &hit)
{
    return Vec3D();
}

RGBColor AmbientLight::incidentRadiance(HitInfo &hit)
{
    return this->lightColor;
}

Point3D AmbientLight::getPos()
{
    return Point3D();
}

RGBColor AmbientLight::getColor()
{
    return this->lightColor;
}

#endif