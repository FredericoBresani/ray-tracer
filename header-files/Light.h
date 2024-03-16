#ifndef __LIGHT__
#define __LIGHT__

#include "Vectors.h"
#include "RGBColor.h"
#include "HitInfo.h"
#include "Points.h"
#include <vector>
#include "GeometricObject.h"

class Light
{
    public:
        virtual Vec3D getDirection(HitInfo &hit) = 0;
        virtual RGBColor incidentRadiance(HitInfo &hit) = 0;
        virtual Point3D getPos() = 0;
        virtual RGBColor getColor() = 0;
        virtual bool castShadows() = 0;
        virtual void sampleLight() = 0;
        virtual bool isExtense() = 0;
        virtual Object* getLightModel() = 0;
        virtual std::vector<Point3D> getLightSamples() = 0;
    protected:
        bool shadows;
        std::vector<Point3D> light_samples;
        int n_samples;
};

#endif