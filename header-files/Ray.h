#ifndef __RAY__
#define __RAY__

#include "Vectors.h"
#include "Points.h"

class Ray {
    public:
        Point3D origin;
        Vec3D direction;
        Ray(void);                              // default constructor
        Ray(const Point3D& o, const Vec3D& dir): origin(o), direction(dir) {}// constructor
        Ray(const Ray& ray);                    // copy constructor
        Ray& operator= (const Ray& rhs);        // assignment operator
        ~Ray() {}                               // destructor
};

#endif