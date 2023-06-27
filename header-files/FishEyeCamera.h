#ifndef __FISHEYE__
#define __FISHEYE__

#include <vector>
#include <fstream>
#include "Camera.h"
#include "GeometricObject.h"
#include "Ambient.h"
#include "Light.h"
#include "Points.h"
#include "Vectors.h"
#include "Sampler.h"
#include "RGBColor.h"
#include "World.h"
#include "Definitions.h"


class FishEyeCamera : public Camera {
    public:
        FishEyeCamera(int _hr, int _vr, double d, const Vec3D& _up, const Point3D& pos, const Point3D &_lookAt, float p, int s, double fAngle)
        {
            hr = _hr;
            vr = _vr;
            distance = d;
            up = _up;
            cameraPos = pos;
            lookAt = _lookAt;
            pixelSize = p;
            samples = s;
            fishAngle = fAngle*M_PI;
        }
        ~FishEyeCamera() {}
    private:
        void set_sampler();
        void render(std::vector<Object*> objetos, std::vector<Light*>& lights, Ambient& ambient);
};

void FishEyeCamera::render(std::vector<Object*> objetos, std::vector<Light*>& lights, Ambient& ambient)
{
    Vec3D toPixel = w*distance + right*(-pixelQtnH/2.0) + iup*(pixelQtnV/2.0);
    Vec3D down;
    Vec3D dir = toPixel;
    float wAngle = fishAngle/2.0;
    float uAngle = fishAngle/2.0;
    std::vector<RGBColor> pixels;
    for (int i = 0; i < pixelQtnH*pixelQtnV; i++) 
    {
        if ((i) % (int)pixelQtnH == 0 && i != 0)
        {
            down = down - iup;
            dir = toPixel + down;
        } else {
            dir = dir + right;
        }
        Point3D screenP = cameraPos + dir;
        Point2D screenCoordinates = this->worldToScreenCoordinates(screenP, cameraPos);
        Vec2D diskVec = Vec2D(screenCoordinates.x, screenCoordinates.y);
        RGBColor sum;
        
        for (int iSamples = 0; iSamples < sampler_ptr->get_num_samples(); iSamples++)
        {
            Point2D aliasUnit = sampler_ptr->sample_unit_square();
            Vec3D sampleX = right*aliasUnit.x;
            Vec3D sampleY = iup*(-1)*(aliasUnit.y);
            Vec3D aliasDir = dir + sampleX + sampleY;
            Point3D aliasPoint = cameraPos + aliasDir;
            Point2D aliasScreenCoordinates = this->worldToScreenCoordinates(aliasPoint, cameraPos);
            Vec2D aliasVec = Vec2D(aliasScreenCoordinates.x, aliasScreenCoordinates.y);
            if (Vec2D::norma(aliasVec) <= 1) {
                float aliasDistance = Vec2D::norma(aliasVec);
                float auxWAngle = wAngle*aliasDistance;
                float cosW = std::cos(auxWAngle); 
                float senW = std::sin(auxWAngle);
                Vec3D direction = u*(aliasScreenCoordinates.x*senW) + v*(senW*aliasScreenCoordinates.y) + w*cosW;
                sum = sum + trace(Ray(cameraPos, direction), objetos, (*this), lights, &ambient, ambient.depth);
            } else {

            }
        }
        pixels.push_back(sum/(double)(sampler_ptr->get_num_samples()));
    }
    std::ofstream pixelOutput("./image.ppm", std::ios::out | std::ios::binary);
    pixelOutput << "P6\n" << pixelQtnH << " " << pixelQtnV << "\n255\n";
    for (int i = 0; i < pixelQtnH*pixelQtnV; i++)
    {
        pixelOutput <<(unsigned char)(std::max(double(1), pixels[i].r)) <<
            (unsigned char)(std::max(double(1), pixels[i].g)) <<
            (unsigned char)(std::max(double(1), pixels[i].b));
    }
    pixelOutput.close();

}

void FishEyeCamera::set_sampler()
{
    sampler_ptr = new JitteredSampler(samples);
}


#endif