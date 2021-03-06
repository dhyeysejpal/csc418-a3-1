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

    Vector3D dir = worldToModel * ray.dir;
    Point3D origin = worldToModel * (ray.origin - ray.time * get_velocity());
    // Define the normal as (0, 0, 1)
    Vector3D normal(0, 0, 1);

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

    Vector3D dir = worldToModel * ray.dir;
    // Randomly skew the origin to simulate the object in motion.
    Point3D origin = worldToModel * (ray.origin - ray.time * get_velocity());

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

        if (t <= 0 || (!ray.intersection.none && t > ray.intersection.t_value)) {
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

bool DysonSphere::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
        const Matrix4x4& modelToWorld ) {

    Vector3D dir = worldToModel * ray.dir;
    Point3D origin = worldToModel * ray.origin;

    if (origin[0] * origin[0] + origin[1] * origin[1] + origin[2] * origin[2] >= 1) {
        // Outside of the sphere, so won't be able to see anything.
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

    // There are always two solutions. Take the one in front.
    t = (-1 * b + sqrt(discriminant)) / (2 * a); // Always positive since we're inside the sphere already.

    if (!ray.intersection.none && t > ray.intersection.t_value) {
        // Intersection point is behind another object.
        return false;
    }

    x = t * dir[0] + origin[0];
    y = t * dir[1] + origin[1];
    z = t * dir[2] + origin[2];

    // The normal points to the origin of the sphere
    Vector3D normal(-x, -y, -z);
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

bool RightCylinder::intersect( Ray3D& ray, const Matrix4x4& worldToModel,
        const Matrix4x4& modelToWorld ) {
    Vector3D dir = worldToModel * ray.dir;
    Point3D origin = worldToModel * (ray.origin - ray.time * get_velocity());
    // Define the normal as (0, 0, 1)
    Vector3D normal;
    Point3D intersection;
    bool res = false;

    // Solve for collisions with the cylinder defined by height |z'| <= 1, radius x**2 + y**2 = 1
    double x, y, z, t;
    double top_t, bottom_t, side_t;

    if (dir[2] != 0) {
        z = 1.0;
        top_t = (z - origin[2]) / dir[2];

        // Calculate point of intersection.
        x = top_t * dir[0] + origin[0]; // x and y are in model coordinates.
        y = top_t * dir[1] + origin[1];

        if (top_t > 0 && pow(x, 2) + pow(y, 2) <= 1) {
            res = true;
            t = top_t;
            normal = Vector3D(0, 0, 1);
            intersection = Point3D(x, y, z);
            ray.intersection.none = false;
        }

        z = -1.0;
        bottom_t = (z - origin[2]) / dir[2];

        // Calculate point of intersection.
        x = bottom_t * dir[0] + origin[0]; // x and y are in model coordinates.
        y = bottom_t * dir[1] + origin[1];

        if (bottom_t > 0 && (!res || bottom_t < t) && pow(x, 2) + pow(y, 2) <= 1) {
            res = true;
            t = bottom_t;
            normal = Vector3D(0, 1, 0);
            intersection = Point3D(x, y, z);
            ray.intersection.none = false;
        }
    }

    // Detect collision with cylinder side.
    // Solve x**2 + y**2 - r**2 = 0
    // = (origin[0] + t*dir[0])**2 + (origin[1] + t*dir[1])**2 - 1
    // = t**2(dir[0]**2 + dir[1]**2) + 2t(origin[0]dir[0] + origin[1]dir[1]) + origin[0]**2 + origin[1]**2 - 1
    double a = pow(dir[0], 2) + pow(dir[1], 2);
    double b = 2 * (dir[0] * origin[0] + dir[1] * origin[1]);
    double c = pow(origin[0], 2) + pow(origin[1], 2) - 1;

    double discriminant = pow(b, 2) -  4 * a * c;

    if (discriminant >= 0) {
        // Intersection point exists. Only look at the closest one, because it's in front and if it's negative,
        // we're inside the cylinder anyway
        side_t = (-b - sqrt(discriminant)) / (2 * a);

        if (side_t > 0 && (!res || side_t )) {
            x = side_t * dir[0] + origin[0]; // x and y are in model coordinates.
            y = side_t * dir[1] + origin[1];
            z = side_t * dir[2] + origin[2];

            if (std::abs(z) <= 1) {
                res = true;
                t = side_t;
                normal = Vector3D(x, y, 0);
                intersection = Point3D(x, y, z);
                ray.intersection.none = false;
            }
        }
    }

    if (res) {
        normal = transNorm(worldToModel, normal);
        normal.normalize();

        ray.intersection.normal = normal;

        // Calculate the intersection point in world coordinates.
        Point3D intersection(x, y, 0);
        intersection = modelToWorld * intersection;

        ray.intersection.point = intersection;

        ray.intersection.t_value = t;

        
    }
    
    

    // Ray doesn't intersect within the unit square's bounds.
    return res;
}