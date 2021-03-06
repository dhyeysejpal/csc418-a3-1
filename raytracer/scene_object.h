/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
            CSC418, SPRING 2005

        classes defining primitives in the scene

***********************************************************/

#include "util.h"

// All primitives should provide a intersection function.  
// To create more primitives, inherit from SceneObject.
// Namely, you can create, Sphere, Cylinder, etc... classes
// here.
class SceneObject {
public:
    // Returns true if an intersection occured, false otherwise.
    virtual bool intersect( Ray3D&, const Matrix4x4&, const Matrix4x4& ) = 0;
    bool in_motion() { return _velocity.length() > 0; }
    void set_velocity(Vector3D val) { _velocity = val; }
    Vector3D get_velocity() { return _velocity; }
private:
    // Instantaneous velocity on the interval [0, t] is what determines motion blur.
    Vector3D _velocity;
};

// Example primitive you can create, this is a unit square on 
// the xy-plane.
class UnitSquare : public SceneObject {
public:
    bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
            const Matrix4x4& modelToWorld );
};

class UnitSphere : public SceneObject {
public:
    bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
            const Matrix4x4& modelToWorld );
};

class DysonSphere : public SceneObject {
public:
    bool intersect( Ray3D& ray, const Matrix4x4& worldToModel,
            const Matrix4x4& modelToWorld );
};

class RightCylinder : public SceneObject {
public:
    bool intersect(Ray3D& ray, const Matrix4x4& worldToModel,
            const Matrix4x4& modelToWorld);
};
