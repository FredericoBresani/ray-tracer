#ifndef __WORLD__
#define __WORLD__


#include <vector>
#include "Definitions.h"
#include "RGBColor.h"
#include "Points.h"
#include "Vectors.h"
#include "HitInfo.h"
#include "Ray.h"
#include "Camera.h"
#include "GeometricObject.h"
#include "Ambient.h"
#include "Light.h"


RGBColor trace(const Point3D& origin, const Point3D& pixel, std::vector<Object*>& objetos, Camera &camera, std::vector<Light*> lights, Ambient *ambient, int depth)
{
    double t = infinity;
    double tmin = infinity;
    double kd, ks, ka, kr, kt, phongExp; 
    HitInfo *hInfo = new HitInfo();
    Ray *ray = new Ray(origin, pixel - origin);
    RGBColor color;
    if (depth == 0) {
        return color;
    }
    for (int i = 0; i < objetos.size(); i++)
    {
        if (objetos[i]->rayObjectIntersect(*ray, &t, *hInfo))
        {
            if (t < tmin)
            {
                tmin = t;
                hInfo->hit_location = ray->origin + ray->direction*tmin;
                hInfo->normal = objetos[i]->getNormal(hInfo->hit_location, *ray);
                color = objetos[i]->getColor();
                kd = objetos[i]->getKd();
                ks = objetos[i]->getKs();
                ka = objetos[i]->getKa();
                kr = objetos[i]->getKr();
                kt = objetos[i]->getKt();
                phongExp = objetos[i]->getPhongExp();
                // color = setPixelColorNormal(hit.normal);
                // color = setPixelColorCoordinates(hit.hit_location);
            }
        }
    }
    
    if (tmin == infinity)
    {
        // return setBackgroundSmoothness(pixel, &camera);
        // return setBackgroundRGBCoordinates(pixel, &camera);
        // return RGBColor(0.0, 0.0, 0.0);
        return RGBColor(255.0, 255.0, 255.0);
        // return RGBColor(190.0, 230.0, 255.0);
    } else {
        double difuseIndice = 0, rMax = 0, gMax = 0, bMax = 0, reflectiveness = 0;
        Vec3D resultingColor, mixedColor;
        RGBColor colorFilter = ambient->color*ka;
        hInfo->toCamera = camera.getPos() - hInfo->hit_location;
        hInfo->toCamera = hInfo->toCamera.normalize(hInfo->toCamera);
        hInfo->viewerReflex = ((hInfo->normal*2)*(hInfo->normal*hInfo->toCamera)) - hInfo->toCamera;
        hInfo->viewerReflex = hInfo->viewerReflex.normalize(hInfo->viewerReflex);
        Point3D hitPoint = hInfo->hit_location + hInfo->normal*0.001;

        Vec3D auxVec = hInfo->viewerReflex^hInfo->normal;

        Vec3D auxReflex = auxVec^hInfo->viewerReflex;
        Point3D reflexPoint = hitPoint + (hInfo->viewerReflex*10.0) + (auxVec*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr)) + ((auxVec*(-1.0))*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr)) + (auxReflex*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr)) + ((auxReflex*(-1.0))*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr));


        Vec3D auxNormal = auxVec^hInfo->normal;
        Point3D difusePoint = hitPoint + hInfo->normal + (auxVec*((double)std::rand()/(double)RAND_MAX)*kd) + ((auxVec*(-1.0))*((double)std::rand()/(double)RAND_MAX)*kd) + (auxNormal*((double)std::rand()/(double)RAND_MAX)*kd) + ((auxNormal*(-1.0))*((double)std::rand()/(double)RAND_MAX)*kd);
        
        reflectiveness = 1.0 + ambient->ir;
        for (int l = 0; l < lights.size(); l++) {
            hInfo->toLight = lights[l]->lightPos - hInfo->hit_location;
            hInfo->toLight = hInfo->toLight.normalize(hInfo->toLight);
            hInfo->reflection = ((hInfo->normal*2)*(hInfo->normal*hInfo->toLight)) - hInfo->toLight;
            hInfo->reflection = hInfo->reflection.normalize(hInfo->reflection);
            mixedColor = Vec3D(lights[l]->lightColor.x*(color.x*reflectiveness), lights[l]->lightColor.y*(color.y*reflectiveness), lights[l]->lightColor.z*(color.z*reflectiveness))/255.0;
            resultingColor = resultingColor + mixedColor*kd*std::max(hInfo->normal*hInfo->toLight, 0.0);
        }
         
        color = RGBColor(colorFilter.x*resultingColor.x, colorFilter.y*resultingColor.y, colorFilter.z*resultingColor.z)/255.0;
        color = RGBColor(std::min(color.x , 255.0), std::min(color.y, 255.0), std::min(color.z, 255.0));
        if (kr > 0) {
            color = color + trace(hitPoint, reflexPoint, objetos, camera, lights, ambient, depth - 1);
            color = color/2.0;
        } else {
            color = color + trace(hitPoint, difusePoint, objetos, camera, lights, ambient, depth - 1);
            color = color/2.0;
        }
    }
    return color;
}


#endif