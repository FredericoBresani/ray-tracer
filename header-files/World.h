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


bool inShadow(Ray ray, std::vector<Object*> &objetos, float lightDistance)
{
    double t = infinity;
    HitInfo *shadowHit = new HitInfo();
    for (int i = 0; i < objetos.size(); i++)
    {
        if (objetos[i]->rayObjectIntersect(ray, &t, *shadowHit)) {
            if (t < lightDistance) {
                return true;
            }
        }
    }
    return false;
}

RGBColor trace(const Ray &ray, std::vector<Object*>& objetos, Camera &camera, std::vector<Light*> lights, Ambient *ambient, int depth)
{
    double t = infinity;
    double tmin = infinity;
    double kd, ks, ka, kr, kt, phongExp; 
    HitInfo *hInfo = new HitInfo();
    RGBColor color, objectColor, flatColor, difuseColor, specularColor, reflectiveColor, transparentColor;
    if (depth == 0) {
        return color;
    }
    for (int i = 0; i < objetos.size(); i++)
    {
        if (objetos[i]->rayObjectIntersect(ray, &t, *hInfo))
        {
            if (t < tmin)
            {
                tmin = t;
                hInfo->hit_location = ray.origin + ray.direction*tmin;
                hInfo->normal = objetos[i]->getNormal(hInfo->hit_location, ray);
                objectColor = objetos[i]->getColor();
                kd = objetos[i]->getKd();
                ks = objetos[i]->getKs();
                ka = objetos[i]->getKa();
                kr = objetos[i]->getKr();
                kt = objetos[i]->getKt();
                phongExp = objetos[i]->getPhongExp();
                hInfo->material_pointer = objetos[i]->material;
                // color = setPixelColorNormal(hit.normal);
                // color = setPixelColorCoordinates(hit.hit_location);
            }
        }
    }
    if (hInfo->hit_object) {
        //shade hit location
        double difuseIndice = 0, rMax = 0, gMax = 0, bMax = 0, reflectiveness = 0;
        RGBColor resultingColor, mixedColor, specularColor;
        hInfo->toCamera = Vec3D::normalize(camera.getPos() - hInfo->hit_location);
        hInfo->viewerReflex = Vec3D::normalize(((hInfo->normal*2)*(hInfo->normal*hInfo->toCamera)) - hInfo->toCamera);
        Point3D hitPoint = hInfo->hit_location + hInfo->normal*0.001;

        Vec3D auxVec = hInfo->viewerReflex^hInfo->normal;

        Vec3D auxReflex = auxVec^hInfo->viewerReflex;
        Vec3D reflexDirection = (hInfo->viewerReflex*10.0) + (auxVec*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr)) + ((auxVec*(-1.0))*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr)) + (auxReflex*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr)) + ((auxReflex*(-1.0))*((double)std::rand()/(double)RAND_MAX)*(1.0 - kr));


        Vec3D auxNormal = auxVec^hInfo->normal;
        Vec3D difuseDirection = hInfo->normal + (auxVec*((double)std::rand()/(double)RAND_MAX)*kd) + ((auxVec*(-1.0))*((double)std::rand()/(double)RAND_MAX)*kd) + (auxNormal*((double)std::rand()/(double)RAND_MAX)*kd) + ((auxNormal*(-1.0))*((double)std::rand()/(double)RAND_MAX)*kd);
        
        reflectiveness = 1.0 + ambient->ir;
        for (int l = 0; l < lights.size(); l++) {
            hInfo->toLight = lights[l]->getDirection((*hInfo));
            // implement the attenuation of light by the distance
            double lightDistance = Vec3D::norma(lights[l]->getPos() - hInfo->hit_location); 
            hInfo->reflection = Vec3D::normalize(((hInfo->normal*2)*(hInfo->normal*hInfo->toLight)) - hInfo->toLight);
            // add some flag to deactiviate shadows
            if (!inShadow(Ray(hInfo->hit_location, hInfo->toLight), objetos, lightDistance))
            {
                mixedColor = ((lights[l]->getColor()^objectColor)*reflectiveness)/255.0;
                specularColor = lights[l]->getColor()*(ks*pow(std::max(hInfo->reflection*hInfo->toCamera, 0.0), phongExp));
                resultingColor = resultingColor + mixedColor*std::max(hInfo->normal*hInfo->toLight, 0.0) + specularColor;
            } 
        }

        // float maxComponent = std::max(resultingColor.r, std::max(resultingColor.g, resultingColor.b));
        flatColor = ambient->color*ka + resultingColor;
        /* if (maxComponent > 255.0) {
            flatColor = (resultingColor/maxComponent)*255.0;
        } else {
            flatColor = resultingColor;
        }*/
        flatColor = RGBColor(std::min(flatColor.r, 255.0), std::min(flatColor.g, 255.0), std::min(flatColor.b, 255.0));
        color = flatColor;

        if (kd > 0) {
            // color = (color + trace(Ray(hitPoint, difuseDirection), objetos, camera, lights, ambient, depth - 1))*(kd/2.0);
        }
        if (kr > 0) {
            color = (color + trace(Ray(hitPoint, reflexDirection), objetos, camera, lights, ambient, depth - 1))*(kr/2.0);
        }

        return color;
    } else {
        //return background color
        // return setBackgroundSmoothness(pixel, &camera);
        // return setBackgroundRGBCoordinates(pixel, &camera);
        // return RGBColor(0.0, 0.0, 0.0);
        return RGBColor(255.0, 255.0, 255.0);
        // return RGBColor(190.0, 230.0, 255.0);
    }
    
}


#endif