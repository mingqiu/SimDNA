#ifndef _VECTOR3D_H
#define _VECTOR3D_H
#include <cmath>
#include <iosfwd>

template<typename T>
class Vector3D {
    T _x, _y, _z;
public:
    Vector3D() : _x{}, _y{}, _z{} {}
    Vector3D( T x_, T y_, T z_ ) : _x{x_}, _y{y_}, _z{z_} {}
    inline T x() const { return _x; }
    inline T y() const { return _y; }
    inline T z() const { return _z; }
    T norm() const;
    T normsq() const;
    friend std::istream& operator>>( std::istream &is, Vector3D<T> &vec ) {
        is >> vec._x >> vec._y >> vec._z;
        return is;
    }
};

template<typename T>
inline T Vector3D<T>::normsq() const {
    return _x * _x + _y * _y + _z * _z;
}

template<typename T>
inline T Vector3D<T>::norm() const {
    return sqrt( normsq() );
}


template<>
inline float Vector3D<float>::norm() const {
    return sqrtf( normsq() );
}

template<typename T>
inline const Vector3D<T> operator+( const Vector3D<T> &a, const Vector3D<T> &b ) {
    return Vector3D<T>{ a.x() + b.x(), a.y() + b.y(), a.z() + b.z() };
}

template<typename T>
inline bool operator==(const Vector3D<T> &thisone, const Vector3D<T> &other)
{ return (
            thisone.x() == other.x()
            && thisone.y() == other.y()
            && thisone.z() == other.z()
    );
}



template<typename T>
inline const Vector3D<T> operator-( const Vector3D<T> &a, const Vector3D<T> &b ) {
    return Vector3D<T>{ a.x() - b.x(), a.y() - b.y(), a.z() - b.z() };
}

// Vector * scalar
template<typename T>
inline const Vector3D<T> operator*( const Vector3D<T> &a, T b ) {
    return Vector3D<T>{ a.x() * b, a.y() * b, a.z() * b };
}

// scalar * Vector
template<typename T>
inline const Vector3D<T> operator*( T a, const Vector3D<T> &b ) {
    return Vector3D<T>{ a * b.x(), a * b.y(), a * b.z() };
}

// Vector / scalar
template<typename T>
inline const Vector3D<T> operator/( const Vector3D<T> &a, T b ) {
    return Vector3D<T>{ a.x() / b, a.y() / b, a.z() / b };
}

// Other useful math
template<typename T>
inline T dot( const Vector3D<T> &a, const Vector3D<T> &b ) {
    return a.x() * b.x() + a.y() * b.y() + a.z() * b.z();
}

template<typename T>
inline Vector3D<T> cross( const Vector3D<T> &a, const Vector3D<T> &b ) {
    return Vector3D<T>{ a.y() * b.z() - a.z() * b.y(),
                        a.z() * b.x() - a.x() * b.z(),
                        a.x() * b.y() - a.y() * b.x() };
}

template<typename T>
inline T cube( T x ) {
    return x * x * x;
}

template<typename T>
std::ostream& operator<<( std::ostream &os, const Vector3D<T> &vec ) {
    os << vec.x() << "\t" << vec.y() << "\t" << vec.z();
    return os;
}

template<typename T>
inline T distSQ(const Vector3D<T> &a, const Vector3D<T> &b) {
    return (a.x() - b.x()) * (a.x() - b.x()) + (a.y() - b.y()) * (a.y() - b.y()) + (a.z() - b.z()) * (a.z() - b.z());
}

template<typename T>
inline T dist(const Vector3D<T> &a, const Vector3D<T> &b) {
    return sqrt(distSQ(a, b));
}

typedef Vector3D<float> Vector3Df;
typedef Vector3D<double> Vector3Dd;
typedef Vector3D<int> Vector3Di;

#endif // _VECTOR3D_H
