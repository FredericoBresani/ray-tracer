#ifndef __POINTLIGHT__
#define __POINTLIGHT__

#include "Light.h"
#include "Points.h"
#include "RGBColor.h"


class PointLight: public Light {
    public:
        Point3D lightPos;
        RGBColor lightColor;
        PointLight(const Point3D& pos, const RGBColor& color): lightPos(pos), lightColor(color) {}
        ~PointLight() {}

        Vec3D getDirection();
        RGBColor incidentRadiance();
};

Vec3D PointLight::getDirection()
{
    return Vec3D();
}

RGBColor PointLight::incidentRadiance()
{
    return RGBColor();
}



#endif