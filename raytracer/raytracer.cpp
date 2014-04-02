/***********************************************************
     Starter code for Assignment 3

     This code was originally written by Jack Wang for
            CSC418, SPRING 2005

        Implementations of functions in raytracer.h, 
        and the main function which specifies the 
        scene to be rendered.   

***********************************************************/


#include "raytracer.h"
#include "bmp_io.h"
#include <cmath>
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <unistd.h>

Raytracer::Raytracer() : _lightSource(NULL) {
    _root = new SceneDagNode();
}

Raytracer::~Raytracer() {
    delete _root;
}

SceneDagNode* Raytracer::addObject( SceneDagNode* parent, 
        SceneObject* obj, Material* mat ) {
    SceneDagNode* node = new SceneDagNode( obj, mat );
    node->parent = parent;
    node->next = NULL;
    node->child = NULL;
    
    // Add the object to the parent's child list, this means
    // whatever transformation applied to the parent will also
    // be applied to the child.
    if (parent->child == NULL) {
        parent->child = node;
    }
    else {
        parent = parent->child;
        while (parent->next != NULL) {
            parent = parent->next;
        }
        parent->next = node;
    }
    
    return node;;
}

LightListNode* Raytracer::addLightSource( LightSource* light ) {
    LightListNode* tmp = _lightSource;
    _lightSource = new LightListNode( light, tmp );
    return _lightSource;
}

void Raytracer::rotate( SceneDagNode* node, char axis, double angle ) {
    Matrix4x4 rotation;
    double toRadian = 2*M_PI/360.0;
    int i;
    
    for (i = 0; i < 2; i++) {
        switch(axis) {
            case 'x':
                rotation[0][0] = 1;
                rotation[1][1] = cos(angle*toRadian);
                rotation[1][2] = -sin(angle*toRadian);
                rotation[2][1] = sin(angle*toRadian);
                rotation[2][2] = cos(angle*toRadian);
                rotation[3][3] = 1;
            break;
            case 'y':
                rotation[0][0] = cos(angle*toRadian);
                rotation[0][2] = sin(angle*toRadian);
                rotation[1][1] = 1;
                rotation[2][0] = -sin(angle*toRadian);
                rotation[2][2] = cos(angle*toRadian);
                rotation[3][3] = 1;
            break;
            case 'z':
                rotation[0][0] = cos(angle*toRadian);
                rotation[0][1] = -sin(angle*toRadian);
                rotation[1][0] = sin(angle*toRadian);
                rotation[1][1] = cos(angle*toRadian);
                rotation[2][2] = 1;
                rotation[3][3] = 1;
            break;
        }
        if (i == 0) {
            node->trans = node->trans*rotation;     
            angle = -angle;
        } 
        else {
            node->invtrans = rotation*node->invtrans; 
        }   
    }
}

void Raytracer::translate( SceneDagNode* node, Vector3D trans ) {
    Matrix4x4 translation;
    
    translation[0][3] = trans[0];
    translation[1][3] = trans[1];
    translation[2][3] = trans[2];
    node->trans = node->trans*translation;  
    translation[0][3] = -trans[0];
    translation[1][3] = -trans[1];
    translation[2][3] = -trans[2];
    node->invtrans = translation*node->invtrans; 
}

void Raytracer::scale( SceneDagNode* node, Point3D origin, double factor[3] ) {
    Matrix4x4 scale;
    
    scale[0][0] = factor[0];
    scale[0][3] = origin[0] - factor[0] * origin[0];
    scale[1][1] = factor[1];
    scale[1][3] = origin[1] - factor[1] * origin[1];
    scale[2][2] = factor[2];
    scale[2][3] = origin[2] - factor[2] * origin[2];
    node->trans = node->trans*scale;    
    scale[0][0] = 1/factor[0];
    scale[0][3] = origin[0] - 1/factor[0] * origin[0];
    scale[1][1] = 1/factor[1];
    scale[1][3] = origin[1] - 1/factor[1] * origin[1];
    scale[2][2] = 1/factor[2];
    scale[2][3] = origin[2] - 1/factor[2] * origin[2];
    node->invtrans = scale*node->invtrans; 
}

Matrix4x4 Raytracer::initInvViewMatrix( Point3D eye, Vector3D view, 
        Vector3D up ) {
    Matrix4x4 mat; 
    Vector3D w;
    view.normalize();
    up = up - up.dot(view)*view;
    up.normalize();
    w = view.cross(up);

    mat[0][0] = w[0];
    mat[1][0] = w[1];
    mat[2][0] = w[2];
    mat[0][1] = up[0];
    mat[1][1] = up[1];
    mat[2][1] = up[2];
    mat[0][2] = -view[0];
    mat[1][2] = -view[1];
    mat[2][2] = -view[2];
    mat[0][3] = eye[0];
    mat[1][3] = eye[1];
    mat[2][3] = eye[2];

    return mat; 
}

void Raytracer::traverseScene( SceneDagNode* node, Ray3D& ray ) {
    SceneDagNode *childPtr;

    // Applies transformation of the current node to the global
    // transformation matrices.
    _modelToWorld = _modelToWorld*node->trans;
    _worldToModel = node->invtrans*_worldToModel; 
    if (node->obj) {
        // Perform intersection.
        if (node->obj->intersect(ray, _worldToModel, _modelToWorld)) {
            ray.intersection.mat = node->mat;
        }
    }
    // Traverse the children.
    childPtr = node->child;
    while (childPtr != NULL) {
        traverseScene(childPtr, ray);
        childPtr = childPtr->next;
    }

    // Removes transformation of the current node from the global
    // transformation matrices.
    _worldToModel = node->trans*_worldToModel;
    _modelToWorld = _modelToWorld*node->invtrans;
}

void Raytracer::computeShading( Ray3D& ray ) {
    Colour col(0.0, 0.0, 0.0);
    for (LightListNode *curLight = _lightSource; curLight != NULL; curLight = curLight->next) {
        if (curLight == NULL) break;
        // Each lightSource provides its own shading function.
        Colour previous_colour(ray.col);
        int i;
        // TODO: Implement shadows here if needed.
        for (i = 0; i < _shadow_samples; i++) {
            // Cast a ray between the intersection point and the light source.
            // If it hits something, the point is in shadow, so shade only
            // using ambient lighting.
            Vector3D l = curLight->light->get_position() - ray.intersection.point;
            double light_dist = l.length();
            l.normalize();
            Ray3D shadow(ray.intersection.point + 0.001 * l, l);
            shadow.num_bounces = ray.num_bounces;
            traverseScene(_root, shadow);

            ray.in_shadow = _shadows_enabled
                && !shadow.intersection.none
                && (shadow.intersection.point - ray.intersection.point).length() < light_dist;
            
            curLight->light->shade(ray);

            if (!_shadows_enabled || !curLight->light->supports_soft_shadows()) {
                // For efficiency, only perform multiple samples for light sources that support it.
                break;
            }
        }

        if (i > 0) {
            // Take the average of the number of samples taken.
            ray.col = (1.0 / i) * ray.col; 
        }
        // TODO: fix shadow balance so that having more samples taken from an AreaLight doesn't
        // cause darker shadows than when fewer samples are taken, under multiple light sources.
        // Colour is simply additive between lights.
        ray.col = ray.col + previous_colour;

        ray.col.clamp();
    }
}

void Raytracer::initPixelBuffer() {
    int numbytes = _scrWidth * _scrHeight * sizeof(unsigned char);
    _rbuffer = new unsigned char[numbytes];
    _gbuffer = new unsigned char[numbytes];
    _bbuffer = new unsigned char[numbytes];
    for (int i = 0; i < _scrHeight; i++) {
        for (int j = 0; j < _scrWidth; j++) {
            _rbuffer[i*_scrWidth+j] = 0;
            _gbuffer[i*_scrWidth+j] = 0;
            _bbuffer[i*_scrWidth+j] = 0;
        }
    }
}

void Raytracer::flushPixelBuffer( char *file_name ) {
    bmp_write( file_name, _scrWidth, _scrHeight, _rbuffer, _gbuffer, _bbuffer );
    delete _rbuffer;
    delete _gbuffer;
    delete _bbuffer;
}

Colour Raytracer::shadeRay( Ray3D& ray ) {
    Colour avg_col(0, 0, 0); // Average colour sample, accounting for motion blur.
    if (getBlurSamples() > 0) {
        double t_interval = 1.0 / getBlurSamples();
        for (int i = 0; i < getBlurSamples(); i++) {
            // Take multiple samples from each ray to account for scene objects in motion.
            ray.col = Colour(0, 0, 0); // Need to start fresh every time.
            ray.intersection.none = true;
            Colour col(0.0, 0.0, 0.0);

            // Linearly interpolate time.
            double r = (rand() % 100) / 100.0;
            double time = (t_interval * (i + r));

            ray.time = time;

            traverseScene(_root, ray); 
        
            // Don't bother shading if the ray didn't hit 
            // anything.
            if (!ray.intersection.none) {
                computeShading(ray);
                col = ray.col;
                // For more perfectly reflective surfaces, reflections are less glossy.
                // Define a region around which reflections are randomly cast.
                double a = 1 - ray.intersection.mat->reflectivity;

                if (ray.intersection.mat->reflectivity > 0 && ray.num_bounces < getRecursiveDepth()) {
                    // Fire off more rays originating at the intersection point.
                    Vector3D n(ray.intersection.normal);
                    Vector3D incident = ray.dir - 2 * (n.dot(ray.dir) * n);
                    incident.normalize();
                    // Distribute rays across a patch perpendicular to the incident ray.
                    Vector3D u = incident.unit_normal();
                    Vector3D v = incident.cross(u);
                    
                    Colour reflection_col(0, 0, 0);
                    int i;
                    for (i = 0; i < _reflection_samples; i++) {
                        // If glossy reflections are enabled, generate orthonormal coordinates
                        // and cast off randomly distributed rays from the plane.
                        Ray3D reflection;
                        if (_glossy_reflections) {
                            double s = a * ((rand() % 100) / 100.0 - 0.5);
                            double t = a * ((rand() % 100) / 100.0 - 0.5);
                            Vector3D r = incident + s * u + t * v;
                            r.normalize();
                            reflection = Ray3D(ray.intersection.point + 0.01 * r, r);
                        } else {
                            reflection = Ray3D(ray.intersection.point + 0.01 * incident, incident);
                        }
                        
                        reflection.num_bounces = ray.num_bounces + 1;
                        reflection_col = reflection_col + ray.intersection.mat->reflectivity * ray.intersection.mat->specular * shadeRay(reflection);
                    }

                    reflection_col = (1.0 / i) * reflection_col;
                    col = col + reflection_col;
                    col.clamp();
                }
            }

            avg_col = avg_col + col;
        }
        // Take the average colour.
        avg_col = (1.0 / getBlurSamples()) * avg_col;
    }
    

    //TODO:
    // You'll want to call shadeRay recursively (with a different ray, 
    // of course) here to implement reflection/refraction effects.

    // For refraction, use Snell's Law:
    // new vector d:
    // d ~= (-c2/c1 * r) + (c2/c1 * cos(th1) - cos(th2)) * N
    // sin(th1) / sin(th2) = c1 / c2

    // Note if c2 > c1, ray bends away from the normal, so you might get a non-acute angle.
    // => total internal reflection, ray doesn't go through to the other side of the 
    // material. So check that incoming angle isn't such that sin*(th1) >= c1/c2


    return avg_col; 
}   

void Raytracer::render( int width, int height, Point3D eye, Vector3D view, 
        Vector3D up, double fov, char* fileName ) {
    Matrix4x4 viewToWorld;
    _scrWidth = width;
    _scrHeight = height;
    double factor = (double(height)/2)/tan(fov*M_PI/360.0);

    initPixelBuffer();
    viewToWorld = initInvViewMatrix(eye, view, up);

    // Construct a ray for each pixel.
    for (int i = 0; i < _scrHeight; i++) {
        for (int j = 0; j < _scrWidth; j++) {
            // Sets up ray origin and direction in view space, 
            // image plane is at z = -1.
            Point3D origin(0, 0, 0);
            Point3D imagePlane;
            Colour col;

            if (getAA()) {
                // If an anti-aliasing mode has been specified, perform stratified sampling.
                int n = getAA();
                for (int p = 0; p < n; p++) {
                    for (int q = 0; q < n; q++) {
                        double r = (rand() % 100) / 100.0;
                        double s = (rand() % 100) / 100.0;

                        // Apply jittering to where the ray is fired off.
                        imagePlane[0] = (-double(width)/2 + j + (p + r) / n)/factor;
                        imagePlane[1] = (-double(height)/2 + i + (q + s) / n)/factor;
                        imagePlane[2] = -1;

                        Vector3D dir = viewToWorld * (imagePlane - origin);
                        dir.normalize();

                        Ray3D ray(viewToWorld * origin, dir);

                        col = col + shadeRay(ray);
                    }
                }

                // Calculate the average colour of the pixel.
                col = (1.0 / (n * n)) * col;

            } else {
                // Fire ray off through the centre of the pixel.
                imagePlane[0] = (-double(width)/2 + 0.5 + j)/factor;
                imagePlane[1] = (-double(height)/2 + 0.5 + i)/factor;
                imagePlane[2] = -1;


                Vector3D dir = viewToWorld * (imagePlane - origin);
                dir.normalize();
                
                Ray3D ray(viewToWorld * origin, dir);
                ray.num_bounces = 1; // A ray bounces once off of the first object it hits,


                col = col + shadeRay(ray);
            }

            _rbuffer[i*width+j] = int(col[0]*255);
            _gbuffer[i*width+j] = int(col[1]*255);
            _bbuffer[i*width+j] = int(col[2]*255);
        }
    }

    flushPixelBuffer(fileName);
}

int main(int argc, char* argv[])
{   
    // Build your scene and setup your camera here, by calling 
    // functions from Raytracer.  The code here sets up an example
    // scene and renders it from two different view points, DO NOT
    // change this if you're just implementing part one of the 
    // assignment.  
    Raytracer raytracer;
    int width = 320; 
    int height = 240;
    raytracer.setAA(0);

    int c;
    extern char *optarg;
    extern int optind, optopt, opterr;

    // Depth of 1 sets rays to emit themselves once and then
    // exit the recursive function.
    raytracer.setRecursiveDepth(0);
    // Raytracer uses at most one shadow sample unless overwritten.
    raytracer.setSS(1);
    raytracer.disableGlossyReflections();
    raytracer.setReflectionSamples(1);
    raytracer.disableShadows();
    raytracer.setBlurSamples(1); // ALways sample an image at at least one point

    while ((c = getopt(argc, argv, ":a:sS:Mw:h:r::g:b:")) != -1) {
        switch (c) {
            case 'M':
                // Auto-medium settings.
                raytracer.setAA(3);
                raytracer.enableShadows();
                raytracer.setSS(10);
                raytracer.setRecursiveDepth(3);
                break;
            case 'a':
                raytracer.setAA(atoi(optarg));
                break;
            case 's':
                raytracer.enableShadows();
                break;
            case 'S':
                // Enable soft shadows.
                raytracer.enableShadows();
                raytracer.setSS(atoi(optarg));
                break;
            case 'r':
                // Enable reflection and refraction.
                if (optarg != NULL) {
                    raytracer.setRecursiveDepth(atoi(optarg));
                } else {
                    // Allow rays to bounce twice by default.
                    raytracer.setRecursiveDepth(3);
                }
                break;
            case 'g':
                // Enable glossy reflections.
                raytracer.enableGlossyReflections();
                raytracer.setReflectionSamples(atoi(optarg));

                break;
            case 'b':
                raytracer.setBlurSamples(atoi(optarg));
                break;
            case 'w':
                width = atoi(optarg);
                break;
            case 'h':
                height = atoi(optarg);
                break;
            default:
                break;
        }
    }

    // Camera parameters.
    Point3D eye(0, 0, 1);
    Vector3D view(0, 0, -1);
    Vector3D up(0, 1, 0);
    double fov = 60;

    // Defines a point light source.
    raytracer.addLightSource( new AreaLight(Point3D(0, 10, 5),
        Vector3D(1, 0, 0), Vector3D(0, 0, 1), Colour(0.9, 0.9, 0.9) ) );
    //raytracer.addLightSource(new PointLight(Point3D(0, 1000, 0), Colour(0.2, 0.2, 0.05)));

    // Defines a material for shading.
    Material gold( Colour(0.24725, 0.1995, 0.0745), Colour(0.75164, 0.60648, 0.22648), 
            Colour(0.628281, 0.555802, 0.366065), 
            51.2 );
    Material jade( Colour(0.135, 0.2225, 0.1575), Colour(0.54, 0.89, 0.63), 
            Colour(0.316228, 0.316228, 0.316228), 
            12.8 );
    Material green_plastic( Colour(0, 0, 0), Colour(0.1, 0.35, 0.1), 
            Colour(0.45, 0.55, 0.45), 
            32, 0.01);
    Material pearl(Colour(0.25, 0.20725, 0.20725), Colour(1.0, 0.829, 0.829),
        Colour(0.296648, 0.296648, 0.296648), 11.264);
     Material mirror(Colour(0, 0, 0), Colour(0, 0, 0),
        Colour(1.0, 1.0, 1.0), 11.264, 1.0);
    Material turquoise(Colour(0.1, 0.18725, 0.1745), Colour(0.396, 0.74151, 0.69102),
        Colour(0.297254, 0.30829, 0.306678), 12.8);
     Material skyblue(Colour(0, 0.74609375, 1), Colour(0, 0, 0), Colour(0, 0, 0), 12.8, 0);

    Material chrome(Colour(0.25, 0.25, 0.25), Colour(0.4, 0.4, 0.4), Colour(0.774597, 0.774597, 0.774597),
        76.8, 0.9);

    // Add a unit square into the scene with material mat.
    // SceneDagNode* sphere = raytracer.addObject( new UnitSphere(), &gold );
    // SceneDagNode* plane = raytracer.addObject( new UnitSquare(), &jade );
    // SceneDagNode *table = raytracer.addObject(new UnitSquare(), &green_plastic);
    
    // Apply some transformations to the unit square.
    double factor1[3] = { 1.0, 2.0, 1.0 };
    double factor2[3] = { 4.0, 4.0, 4.0 };
    double table_factor[3] = { 24.0, 72.0, 1.0 };
    double pool_ball_size[3] = { 3, 3, 3 };
    double sky_size[3] = { 1500, 1500, 1500 };

    // SceneDagNode *sky = raytracer.addObject(new DysonSphere(), &turquoise);
    //     raytracer.scale(sky, Point3D(0, 0, 0), sky_size);

    SceneDagNode *ball1 = raytracer.addObject(new RightCylinder(), &pearl);
        raytracer.translate(ball1, Vector3D(0, -2, 6));
        raytracer.scale(ball1, Point3D(0, 0, 0), pool_ball_size);
    // SceneDagNode *ball2 = raytracer.addObject(new UnitSphere(), &pearl);
    //     raytracer.translate(ball2, Vector3D(1, -2, 10));
    //     raytracer.scale(ball2, Point3D(0, 0, 0), pool_ball_size);
    // SceneDagNode *ball3 = raytracer.addObject(new UnitSphere(), &mirror);
    //     raytracer.translate(ball3, Vector3D(-1, -2, 10));
    //     raytracer.scale(ball3, Point3D(0, 0, 0), pool_ball_size);
    // SceneDagNode *ball4 = raytracer.addObject(new UnitSphere(), &chrome);
    //     raytracer.translate(ball4, Vector3D(2, -2, 14));
    //     raytracer.scale(ball4, Point3D(0, 0, 0), pool_ball_size);
    // SceneDagNode *ball5 = raytracer.addObject(new UnitSphere(), &mirror);
    //     raytracer.translate(ball5, Vector3D(0, -2, 14));
    //     raytracer.scale(ball5, Point3D(0, 0, 0), pool_ball_size);
    // SceneDagNode *ball6 = raytracer.addObject(new UnitSphere(), &pearl);
    //     raytracer.translate(ball6, Vector3D(-2, -2, 14));
    //     raytracer.scale(ball6, Point3D(0, 0, 0), pool_ball_size);
    //     ball6->obj->set_velocity(Vector3D(0, 2.0, 0));
    // SceneDagNode *ball7 = raytracer.addObject(new UnitSphere(), &pearl);
    //     raytracer.translate(ball7, Vector3D(3, -2, 18));
    //     raytracer.scale(ball7, Point3D(0, 0, 0), pool_ball_size);
    // SceneDagNode *ball8 = raytracer.addObject(new UnitSphere(), &pearl);
    //     raytracer.translate(ball8, Vector3D(-3, -2, 18));
    //     raytracer.scale(ball8, Point3D(0, 0, 0), pool_ball_size);
    // SceneDagNode *cue = raytracer.addObject(new UnitSphere(), &pearl);
    //     raytracer.translate(cue, Vector3D(0, -2, 32));
    //     raytracer.scale(cue, Point3D(0, 0, 0), pool_ball_size);

    // SceneDagNode *stick = raytracer.addObject(new RightCylinder(), &gold);

    // raytracer.translate(sphere, Vector3D(0, 0, -5));    
    // raytracer.rotate(sphere, 'x', -45); 
    // raytracer.rotate(sphere, 'z', 45); 
    // raytracer.scale(sphere, Point3D(0, 0, 0), factor1);

    // raytracer.translate(table, Vector3D(0, -3, 0));
    // raytracer.rotate(table, 'x', 270);
    // raytracer.scale(table, Point3D(0, 0, 0), table_factor);

    // raytracer.translate(plane, Vector3D(0, 0, -7)); 
    // raytracer.rotate(plane, 'z', 45);
    // raytracer.scale(plane, Point3D(0, 0, 0), factor2);
    
    // Render it from a different point of view.
    Point3D eye2(20, 10, 28);
    Point3D eye3(2, 0, 18);
    Vector3D view2(-4, -2, -6);

    //raytracer.render(width, height, eye, view, up, fov, (char *) "img1.bmp");
    //raytracer.render(width, height, eye2, view2, up, fov, (char *) "img2.bmp");
    raytracer.render(width, height, eye3, view2, up, fov, (char *) "img2.bmp");
    
    return 0;
}

