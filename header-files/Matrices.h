#ifndef __MATRICES__
#define __MATRICES__

template <typename T> class Matrix4
{  
    public:
        std::vector<std::vector<T>> matrix;
        Matrix4(): matrix(std::vector<std::vector<T>>(4, std::vector<T>(4))) {}
};

template <typename T> class Matrix3
{
    public:
        std::vector<std::vector<T>> matrix;
        Matrix3(): matrix(std::vector<std::vector<T>>(3, std::vector<T>(3))) {}
};

typedef Matrix4<double> Matrix4D;
typedef Matrix3<double> Matrix3D;

#endif