#include <stdio.h>
#include <iostream>
#include <stdlib.h>

// lets use doubles for object-ray intersection and floats for shading calculations

template <typename T> class Vec2
{
    public:
        T x, y;
        Vec2(); // build all the vec2 operations inside this public section
};

template <typename T> class Vec3
{
    public:
        T x, y, z;
        Vec3(); //build all the vec3 operation inside this public section
};

typedef Vec2<float> vec2f;
typedef Vec3<float> vec3f;

class Object 
{
    public:
        Object();
        virtual int rayObjectIntersect();
};

class Sphere: public Object 
{

};

class Triangle: public Object 
{

};

class Plane: public Object 
{

};





int main()
{
    std::cout << "Hello my mans!";
    int test = 0;
    return 0;
}