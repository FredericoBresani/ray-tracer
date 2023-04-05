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
    Vec3D dir, down;
    float uAngle = fishAngle;
    float wAngle = fishAngle/2.0;
    float auxWAngle = fishAngle/pixelQtnV;
    std::vector<RGBColor> pixels;
    for (int i = 0; i < pixelQtnH*pixelQtnV; i++) 
    {
        if ((i) % (int)pixelQtnH == 0)
        {
            down = down - iup;
            dir = toPixel + down;
            if (wAngle < 0)
            {
                wAngle = 2*M_PI - auxWAngle;
                auxWAngle += fishAngle/pixelQtnV;
            } else {
                wAngle -= fishAngle/pixelQtnV;
            }
            uAngle = fishAngle;
        } else {
            dir = dir + right;
            uAngle -= uAngle/pixelQtnH;
        }
        Point3D screenP = cameraPos + dir;
        Point2D screenCoordinates = this->worldToScreenCoordinates(screenP, cameraPos);
        Vec2D diskVec = Vec2D(screenCoordinates.x, screenCoordinates.y);
        Vec3D direction = u*std::cos(uAngle) + v*std::sin(wAngle) + w*std::cos(wAngle) + w*std::sin(uAngle);

        Vec3D sum;
        if (diskVec.norma(diskVec) <= 1) {
            for (int iSamples = 0; iSamples < sampler_ptr->get_num_samples(); iSamples++)
            {
                Point2D aliasUnit = sampler_ptr->sample_unit_square();
                Vec3D sampleX = right*aliasUnit.x;
                Vec3D sampleY = iup*(-1)*(aliasUnit.y);
                Vec3D aliasDir = dir + sampleX + sampleY;
                Point3D aliasPoint = cameraPos + aliasDir;


                sum = sum + trace(Ray(aliasPoint, direction), objetos, (*this), lights, &ambient, ambient.depth);
            }
        } else {
            sum = Vec3D(0.0, 0.0, 0.0);
        }
        pixels.push_back(sum/(double)(sampler_ptr->get_num_samples()));
    }
    std::ofstream pixelOutput("./image.ppm", std::ios::out | std::ios::binary);
    pixelOutput << "P6\n" << pixelQtnH << " " << pixelQtnV << "\n255\n";
    for (int i = 0; i < pixelQtnH*pixelQtnV; i++)
    {
        pixelOutput <<(unsigned char)(std::max(double(1), pixels[i].x)) <<
            (unsigned char)(std::max(double(1), pixels[i].y)) <<
            (unsigned char)(std::max(double(1), pixels[i].z));
    }
    pixelOutput.close();

}

void FishEyeCamera::set_sampler()
{
    sampler_ptr = new JitteredSampler(samples);
}


#endif