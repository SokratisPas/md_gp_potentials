#pragma once

#include <cmath>
#include <stdexcept>

struct Vec3 {
    double x, y, z;

    Vec3() : x(0), y(0), z(0)
    {}
    Vec3(double x, double y, double z) : x(x), y(y), z(z)
    {}
    Vec3(double val) : x(val), y(val), z(val)
    {}

    Vec3 operator+(const Vec3& other) const {
        return Vec3(x + other.x, y + other.y, z + other.z);
    }

    Vec3 operator-(const Vec3& other) const {
        return Vec3(x - other.x, y - other.y, z - other.z);
    }

    Vec3 operator*(double scalar) const {
        return Vec3(x * scalar, y * scalar, z * scalar);
    }

    Vec3 operator/(double scalar) const {
        if (scalar == 0.0) {
            throw std::runtime_error("Division by zero in Vec3 operator/");
        }

        return Vec3(x / scalar, y / scalar, z / scalar);
    }

    double dot(const Vec3& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    double length() const {
        return std::sqrt(dot(*this));
    }
};