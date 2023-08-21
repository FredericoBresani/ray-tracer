#ifndef __LINES__
#define __LINES__

#include "GeometricObject.h"
#include "Vectors.h"
#include "Points.h"
#include "Definitions.h"
#include "Ray.h"
#include "RGBColor.h"
#include "Material.h"

class Line: public Object 
{
    public:
        Point3D origin;
        Vec3D direction, normal;
        RGBColor color;
        double difuseK, specularK, ambientK, reflectionK, transmissionK, phongExp;
        Material *material;
        bool castShadows = false;
        Line(const Point3D &o, const Vec3D &d, const RGBColor &RGB): origin(o), direction(d), color(RGB) {}
        ~Line() {}
        bool rayObjectIntersect(const Ray &ray, double *tmin, HitInfo& hit)
        { 
            double t = ((ray.origin - origin)*normal)/((normal*ray.direction)*(-1.0));
            Point3D location = ray.origin + (ray.direction*t);
            if (location.x > 7 || location.x < -7) return false;
            if (location.y > 7 || location.y < -7) return false;
            if (location.z > 7 || location.z < -7) return false;
            Vec3D toLoc = Vec3D::normalize(location - origin);
            double cos = toLoc*(Vec3D::normalize(direction));
            if (t > kEpsilon && t < (*tmin) && (cos >= -1.0 - 0.000002 && cos <= -1.0 + 0.000002 || cos >= 1.0 - 0.000002 && cos <= 1.0 + 0.000002))
            {
                RGBColor color = this->getColor();
                (*tmin) = t;
                // hit.hit_location = Point3D(color.x, color.y, color.z);
                // hit.normal = Vec3D::normalize(direction);
                return true;
            }
            return false;
        }
        void setNormal(const Point3D& C)
        {
            double closerT = ((origin - C)*direction*(-1.0))/(direction*direction);
            Point3D closerP = origin + (direction*closerT);
            this->normal = C - closerP;
        }
        RGBColor getColor() const
        {
            return this->color;
        }
        double getKd() const
        {
            return this->difuseK;
        }
        double getKs() const
        {
            return this->specularK;
        }
        double getKa() const
        {
            return this->ambientK;
        }
        double getKr() const
        {
            return this->reflectionK;
        }
        double getKt() const
        {
            return this->transmissionK;
        }
        double getPhongExp() const
        {
            return this->phongExp;
        }
        double getIor() const 
        {
            return this->material->ior;
        }
        Vec3D getNormal(const Point3D &hit, const Ray &ray) const
        {
            return Vec3D();
        }
        bool getShadows() const
        {
            return false;
        }
        bool getCastShadows() const
        {
            return false;
        }
};

#endif