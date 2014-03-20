/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
		    CSC418, SPRING 2005

		implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

double max(double a, double b) {
	return (a > b ? a : b);
}

void PointLight::shade(Ray3D& ray) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  

	// Shade the ambient component.
	Colour ambient(ray.intersection.mat->ambient);
	ambient = ambient * _col_ambient;

	// Ambient term shouldn't be additive. Set it to whatever light gives the highest ambient term.
	if (ray.col[0] < ambient[0] && ray.col[1] < ambient[1] && ray.col[2] < ambient[2]) {
		ray.col = ambient;
	}

	// Calculate the diffuse component.
	Vector3D s = (this->get_position() - ray.intersection.point);
	s.normalize();
	Vector3D n = ray.intersection.normal;
	n.normalize();
	Colour diffuse = max(0, s.dot(n)) * ray.intersection.mat->diffuse;
	diffuse = diffuse * _col_diffuse;

	// Calculate the specular component.
	Vector3D b = (ray.origin - ray.intersection.point);
	b.normalize();
	Vector3D r = -1 * s + 2 * (n.dot(s) * n);
	r.normalize();
	Colour specular = pow(max(0, r.dot(b)), ray.intersection.mat->specular_exp) * ray.intersection.mat->specular;
	specular = specular * _col_specular;

	ray.col = ray.col + diffuse + specular;
	ray.col.clamp();
}
