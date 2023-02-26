#ifndef __HITINFO__
#define __HITINFO__

#include "Vectors.h"
#include "Points.h"

class HitInfo {
    public:
        bool hit_object;
        Point3D hit_location;
        Vec3D normal, toLight, toCamera, reflection, refraction, viewerReflex;
        HitInfo() {}
        ~HitInfo() {}
};

#endif