#ifndef __TRIANGLE__
#define __TRIANGLE__

#include "GeometricObject.h"
#include "Points.h"
#include "Vectors.h"
#include "Definitions.h"
#include "Plane.h"
#include "RGBColor.h"
#include "Material.h"

class Triangle: public Object 
{
    public:
        Point3D A;
        Point3D B;
        Point3D C;
        Material *material;
        bool castShadows;
        Triangle(const Point3D &a, const Point3D &b, const Point3D &c, Material *m): A(a), B(b), C(c), material(m) {}
        ~Triangle() {}
        bool rayObjectIntersect(const Ray& ray, double *tmin, HitInfo& info)
        {
            Vec3D tPlaneNormal = (this->A - this->B) ^ (this->A - this->C);
            if (((ray.direction * (-1)) * tPlaneNormal) <= 0)
            {
                tPlaneNormal = (this->A - this->C) ^ (this->A - this->B);
            }
            Material *tempMaterial = new Material{
                color, 0, 0, 0, 0, 0, 0, false, 0
            };
            Plane *tPlane = new Plane(tPlaneNormal, this->A, tempMaterial, false);
            Point3D pHit;
            if (tPlane->rayObjectIntersect(ray, tmin, info)) 
            {   
                free(tPlane);
                if ((*tmin) < kEpsilon) return true;
                pHit = ray.origin + ray.direction*(*tmin);
                //|A.x B.x C.x||a|   |X|
                //|A.y B.y C.y||b| = |Y|
                //|A.z B.z C.z||g|   |Z|
                double det = A.x*B.y*C.z - A.x*C.y*B.z - B.x*A.y*C.z + B.x*C.y*A.z + C.x*A.y*B.z - C.x*B.y*A.z;

                double alpha = (pHit.x*B.y*C.z - pHit.x*C.y*B.z - B.x*pHit.y*C.z + B.x*C.y*pHit.z + C.x*pHit.y*B.z - C.x*B.y*pHit.z)/det;
                if (alpha > 1 || alpha < 0) return false;

                double beta = (A.x*pHit.y*C.z - A.x*C.y*pHit.z - pHit.x*A.y*C.z + pHit.x*C.y*A.z + C.x*A.y*B.z - C.x*pHit.y*A.z)/det;
                if (beta > 1 || beta < 0) return false;

                if (alpha + beta > 1 || alpha + beta < 0) return false;

                // double gama = (A.x*B.y*pHit.z - A.x*pHit.y*B.z - B.x*A.y*pHit.z + B.x*pHit.y*A.z + pHit.x*A.y*B.z - pHit.x*B.y*A.z)/det; 
                
                info.hit_object = true;
                return true;
                
            } else {
                (*tmin) = infinity;
                return false;
            }
        }
        RGBColor getColor() const
        {
            return this->material->color;
        }
        double getKd() const
        {
            return this->material->difuseK;
        }
        double getKs() const
        {
            return this->material->specularK;
        }
        double getKa() const
        {
            return this->material->ambientalK;
        }
        double getKr() const
        {
            return this->material->reflectiveK;
        }
        double getKt() const
        {
            return this->material->transmissionK;
        }
        double getPhongExp() const
        {
            return this->material->roughK;
        }
        double getIor() const 
        {
            return this->material->ior;
        }
        Vec3D getNormal(const Point3D &hit, const Ray &ray) const
        {
            Vec3D tPlaneNormal = (this->A - this->B) ^ (this->A - this->C);
            if (((ray.direction * (-1)) * tPlaneNormal) <= 0)
            {
                tPlaneNormal = (this->A - this->C) ^ (this->A - this->B);
            }
            return Vec3D::normalize(tPlaneNormal);
        }
        bool getShadows() const
        {
            return this->material->getShadows;
        }
        bool getCastShadows() const
        {
            return this->castShadows;
        }
};

#endif