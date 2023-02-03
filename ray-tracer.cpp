#include <stdio.h>
#include <iostream>
#include <vector>

template <typename T>
class Vec2
{
    T x, y;
    Vec2() : x(T(0)), y(T(0)) {}
    Vec2(T _x) : x(T(_x)), y(T(_x)) {}
};

template <typename T> 
class Vec3
{
    T x, y, z;
    Vec3(): x(T(0)), y(T(0)), z(T(0)) {}
    Vec3(T _x): x(T(_x)), y(T(_x)), z(T(_x)) {}
};

typedef Vec2<float> Vec2f;
typedef Vec3<float> Vec3f;


void render()
{
    int hw = 200;
    int vw = 100;
    std::vector<Vec3f> pixelsRGB;
    std::cout << "P3\n" << hw << " " << vw << "\n255\n";
    for (int )
}


int main()
{
    render();
    return 0;
}