/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
            CSC418, SPRING 2005

           light source classes

***********************************************************/

#include "util.h"

// Base class for a light source.  You could define different types
// of lights here, but point light is sufficient for most scenes you
// might want to render.  Different light sources shade the ray 
// differently.
class LightSource {
public:
    virtual void shade( Ray3D& ) = 0;
    virtual Point3D get_position() const = 0;
};

// A point light is defined by its position in world space and its
// colour.
class PointLight : public LightSource {
public:
    PointLight( Point3D pos, Colour col ) : _pos(pos), _col_ambient(col), 
    _col_diffuse(col), _col_specular(col) {}
    PointLight( Point3D pos, Colour ambient, Colour diffuse, Colour specular ) 
    : _pos(pos), _col_ambient(ambient), _col_diffuse(diffuse), 
    _col_specular(specular) {}
    void shade( Ray3D& ray );
    Point3D get_position() const { return _pos; }
    
private:
    Point3D _pos;
    Colour _col_ambient;
    Colour _col_diffuse; 
    Colour _col_specular; 
};

/**
  * Area lighting that is defined by its centre and two vectors representing
  * the patch in space where this light exists.
  */
class AreaLight : public LightSource {
public:
    AreaLight(Point3D centre, Vector3D a, Vector3D b, Colour col)
        : _centre(centre), _a(a), _b(b), _n(a.cross(b)), _col_ambient(col),
            _col_diffuse(col), _col_specular(col) {}
    AreaLight(Point3D centre, Vector3D a, Vector3D b, Colour ambient, Colour diffuse, Colour specular)
        : _centre(centre), _a(a), _b(b), _n(a.cross(b)), _col_ambient(ambient), _col_diffuse(diffuse),
            _col_specular(specular) {}

    void shade(Ray3D &ray);
    Point3D get_position() const { return _centre; }
    Vector3D get_a() const { return _a; }
    Vector3D get_b() const { return _b; }
    Vector3D get_normal() const { return _n; }
private:
    Point3D _centre;
    Vector3D _a;
    Vector3D _b;
    Vector3D _n;
    Colour _col_ambient;
    Colour _col_diffuse;
    Colour _col_specular;
};
