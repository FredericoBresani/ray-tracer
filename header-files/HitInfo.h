#ifndef __HITINFO__
#define __HITINFO__

#include "Vectors.h"
#include "Points.h"

class Material;

class HitInfo {
    public:
        bool hit_object = false;
        Point3D hit_location;
        Vec3D normal, toLight, toCamera, reflection, refraction, viewerReflex;
        Material *material_pointer;
        HitInfo() {}
        ~HitInfo() {}
};

#endif