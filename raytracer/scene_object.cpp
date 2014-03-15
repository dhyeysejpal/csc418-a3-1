/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
            CSC418, SPRING 2005

        implements scene_object.h

***********************************************************/

#include <cmath>
#include <iostream>
#include "scene_object.h"

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
    Point3D eye = worldToModel * ray.origin;

    // Solve: (x', y', z') = eye + t(dir), -0.5 <= x, y <= 0.5, z = 0
    // => x' - eye.x = tx; y' - eye.y = ty; z' -eye.z = tz

    if (dir[2] == 0) {
        // ray in the xy plane, so even if it intersects, the square won't be visible.
        ray.intersection.none = true;
        return false;
    }
    // Solve z' = eye.z + t * z
    double x, y, t;
    t = (0 - eye[2]) / dir[2];
    x = t * dir[0] + eye[0]; 
    y = t * dir[1] + eye[1];

    ray.intersection.none = !(x >= -0.5 && x <= 0.5 && y >= -0.5 && y <= 0.5);

    if (!ray.intersection.none) {
        // Define the normal as (0, 0, 1)
        Vector3D normal(0, 0, 1);
        // Convert it to world space.
        normal = transNorm(worldToModel, normal);
        ray.intersection.normal = normal;

        // Calculate what the t value is for the intersection point in
        // world space.
        Point3D intersection(x, y, 0);
        intersection = modelToWorld * intersection;

        // Pick a component of the vector that is not 0 and solve for t.
        if (ray.dir[2] != 0) {
            t = (intersection[2] - ray.origin[2]) / ray.dir[2];
        } else if (ray.dir[1] != 0) {
            t = (intersection[1] - ray.origin[1]) / ray.dir[1];
        } else {
            // A ray cannot have |v| = 0, so some component must not be 0.
            t = (intersection[0] - ray.origin[0]) / ray.dir[0];
        }

        ray.intersection.t_value = t;
    }

    // Convert ray into object space.

    // Solve for intersection at x = +-0.5, y = +-0.5, z = +-0.5,
    // and pick the intersection point with the lowest t_value.

    // Convert intersection point back into world space.

    return !ray.intersection.none;
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
    
    return false;
}

