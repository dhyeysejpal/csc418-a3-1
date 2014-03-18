/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
            CSC418, SPRING 2005

        implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"
#include <stdio.h>

bool UnitSquare::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
        const Matrix4x4& modelToWorld ) {
    // TODO: implement intersection code for UnitSquare, which is
    // defined on the xy-plane, with vertices (0.5, 0.5, 0), 
    // (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
    // (0, 0, 1).
    //
    // Your goal here is to fill ray.intersection with correct values
    // should an intersection occur.  This includes intersection.point, 
    // intersection.normal, intersection.none, intersection.t_value.   
    //
    // HINT: Remember to first transform the ray into object space  
    // to simplify the intersection test.

    Vector3D dir = worldToModel * ray.dir;
    Point3D origin = worldToModel * ray.origin;

    // Solve: (x', y', z') = origin + t(dir), -0.5 <= x, y <= 0.5, z = 0
    // => x' - origin.x = tx; y' - origin.y = ty; z' -origin.z = tz

    if (dir[2] == 0) {
        // ray is in the xy plane, so even if it intersects, the square won't be visible.
        return false;
    }
    // Solve z' = origin.z + t * z
    double x, y, t;
    t = (0 - origin[2]) / dir[2];

    if (t <= 0 || (!ray.intersection.none && t > ray.intersection.t_value)) {
        // Looking away from the plane or intersection point is behind an
        // earlier intersection, so don't update the intersection.
        return false;
    }

    // Calculate point of intersection.
    x = t * dir[0] + origin[0]; // x and y are in model coordinates.
    y = t * dir[1] + origin[1];

    if (x >= -0.5 && x <= 0.5 && y >= -0.5 && y <= 0.5) {
        // Define the normal as (0, 0, 1)
        Vector3D normal(0, 0, 1);
        // Convert it to world space.
        normal = transNorm(worldToModel, normal);
        normal.normalize();
        ray.intersection.normal = normal;

        // Calculate the intersection point in world coordinates.
        Point3D intersection(x, y, 0);
        intersection = modelToWorld * intersection;

        ray.intersection.point = intersection;

        ray.intersection.t_value = t;

        ray.intersection.none = false;

        return true;
    }

    // Ray doesn't intersect within the unit square's bounds.
    return false;
}

bool UnitSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
        const Matrix4x4& modelToWorld ) {
    // TODO: implement intersection code for UnitSphere, which is centred 
    // on the origin.  
    //
    // Your goal here is to fill ray.intersection with correct values
    // should an intersection occur.  This includes intersection.point, 
    // intersection.normal, intersection.none, intersection.t_value.   
    //
    // HINT: Remember to first transform the ray into object space  
    // to simplify the intersection test.

    Vector3D dir = worldToModel * ray.dir;
    Point3D origin = worldToModel * ray.origin;

    if (origin[0] * origin[0] + origin[1] * origin[1] + origin[2] * origin[2] <= 1) {
        // Inside of the sphere, so won't be able to see anything.
        return false;
    }

    double x, y, z, t;

    Vector3D radial = origin - Point3D(0, 0, 0);
    // Solve for t = 1/dir.dot(dir) * (-dir.dot(radial))
        // +- sqrt(dir.dot(radial)**2 - (dir.dot(dir) * (radial.dot(radial) - r ** 2))
    double a, b, c;
    a = dir.dot(dir);
    b = 2 * dir.dot(radial);
    c = radial.dot(radial) - 1;

    double discriminant = b * b - 4 * a * c;

    if (discriminant >= 0) {
        // There is at least one solution, or at most two.
        // But always take the closest one, 
        t = (-1 * b - sqrt(discriminant)) / (2 * a);

        if (!ray.intersection.none && t > ray.intersection.t_value) {
            // Intersection point is behind another object.
            return false;
        }

        x = t * dir[0] + origin[0];
        y = t * dir[1] + origin[1];
        z = t * dir[2] + origin[2];

        // The normal points straight out from the sphere.
        Vector3D normal(x, y, z);
        // Transform to world space.
        normal = transNorm(worldToModel, normal);
        normal.normalize();
        ray.intersection.normal = normal;

        Point3D intersection(x, y, z);
        ray.intersection.point = modelToWorld * intersection;

        ray.intersection.t_value = t;

        ray.intersection.none = false;

        return true;
    }
    
    // No intersection with the sphere.
    return false;
}
