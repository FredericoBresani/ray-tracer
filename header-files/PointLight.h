#ifndef __POINTLIGHT__
#define __POINTLIGHT__

#include "Light.h"
#include "Points.h"
#include "RGBColor.h"
#include "HitInfo.h"


class PointLight: public Light {
    public:
        PointLight(const Point3D& pos, const RGBColor& color): lightPos(pos), lightColor(color) {}
        ~PointLight() {}

        Vec3D getDirection(HitInfo &hit);
        RGBColor incidentRadiance(HitInfo &hit);
        Point3D getPos();
        RGBColor getColor();
    private:
        Point3D lightPos;
        RGBColor lightColor;
        
};

Vec3D PointLight::getDirection(HitInfo &hit)
{
    Vec3D toLight = this->lightPos - hit.hit_location;
    return toLight.normalize(toLight);
}

RGBColor PointLight::incidentRadiance(HitInfo &hit)
{
    return lightColor;
}

Point3D PointLight::getPos()
{
    return this->lightPos;
}

RGBColor PointLight::getColor()
{
    return this->lightColor;
}



#endif