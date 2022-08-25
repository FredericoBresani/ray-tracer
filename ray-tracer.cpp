#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include <cstdio>

// lets use doubles for object-ray intersection and floats for shading calculations

template <typename T> class Vec2
{
    public:
        T x, y;
        Vec2(): x(T(0)), y(T(0)) {} // build all the vec2 operations inside this public section
        Vec2(T _x): x(T(_x)), y(T(_x)) {}
        Vec2(T _x, T _y): x(T(_x)), y(T(_y)) {}
        // soma e subtração por escalar
        // multiplicacao e divisão por escalar pela esquerda ou direita
        // produto escalar
        // produto vetorial
        // soma vetorial
        // normalização
        // norma
        // multiplicacao por matriz 4x4 a esquerda
        // ponto - ponto
};

template <typename T> class Vec3
{
    public:
        T x, y, z;
        Vec3(): x(T(0)), y(T(0)), z(T(0)) {} //build all the vec3 operation inside this public section
        Vec3(T _x): x(T(_x)), y(T(_x)), z(T(_x)) {}
        Vec3(T _x, T _y, T _z): x(T(_x)), y(T(_y)), z(T(_z)) {}
        // soma e subtração por escalar
        // multiplicacao e divisão por escalar pela esquerda ou direita
        // produto escalar
        // produto vetorial
        // soma vetorial
        // normalização
        // norma
        // multiplicacao por matriz 4x4 a esquerda
        // ponto - ponto
};

template <typename T> class Point2
{
    
};

template <typename T> class Point3
{
    public:
        T x, y, z;
        Point3(): x(T(0)), y(T(0)), z(T(0)) {}
        Point3(T _x): x(T(_x)), y(T(_x)), z(T(_x)) {}
        Point3(T _x, T _y, T _z): x(T(_x)), y(T(_y)), z(T(_z)) {}
        // ponto + vetor, ponto - vetor
        // ponto multiplicado por escalar pela esquerda ou direita
        // soma varicentrica
        // interpolação
        // ponto multiplicado por matriz 4x4 pela esquerda
        // distancia entre dois pontos
};

template <typename T> class Matrix4
{  
    public:
        std::vector<std::vector<T>> matrix;
        Matrix4(): matrix(std::vector<std::vector<T>>(4, std::vector<T>(4))) {}
};

typedef Vec2<double> vec2f;
typedef Vec3<double> vec3f;
typedef Point3<double> point3f;
typedef Matrix4<double> Matrix4f;

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
    Matrix4f matriz = Matrix4f();
    char objectType;
    int _1, _2, _3, _4, _5, _6, _7, _8, _9, _10, _11, _12;
    while (scanf("%c %i %i %i %i %i %i %i %i %i %i %i %i\n", &objectType, &_1, &_2, &_3, &_4, &_5, &_6, &_7, &_8, &_9, &_10, &_11, &_12) != EOF) 
    {
        printf("%c", objectType);
    }
    int test = 1;
    return 0;
}