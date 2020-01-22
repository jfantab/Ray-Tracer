#define GL_SILENCE_DEPRECATION
#include <GLUT/glut.h>
#include <cmath>
#include <iostream>

struct Sphere {
    float center[3];
    float radius;
    float color[3];
    bool isReflexive() {
        return false;
    }

};

struct Triangle {
    float p0[3];
    float p1[3];
    float p2[3];
    float color[3];
    bool isReflexive() {
        return false;
    }

};

struct Ray {
    float p[3];
    float d[3];
};

float pix_col[3] = {0, 0, 0};
Sphere spheres[1];
Triangle triangles[1];
float e1[3], e2[3], e3[3];
float light_pos[3] = {10.0, 30.0, 32.0};
float light_color[3] = {1.0, 1.0, 1.0};

void initSpheres();
void initTriangles();
void initRay();
void drawPixels();
void rayTrace(float[], float[]);
float findSphereIntersection(float[], float[]);
float findTriangleIntersection(float[], float[]);
void calculateLighting(float[], float[], float[], float[]);
float dotProduct(float[], float[]);
void crossProduct(float [], float [], float (&c)[3]);
float findMagnitude(float, float, float);
void normalize(float(&v)[3]);

int main(int argc, char * argv[]) {
    initSpheres();
    initTriangles();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("Ray Tracer");
    glutPositionWindow(500, 300);
    glutReshapeWindow(256, 256);
    glClearColor(0.0, 0.0, 0.0, 1.0);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glOrtho(-128.0, 128.0, -128.0, 128.0, -10.0, 10.0);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glutDisplayFunc(drawPixels);
    glutMainLoop();
    return 0;
}

void initSpheres() {
    spheres[0].center[0] = 0.0;
    spheres[0].center[1] = 50.0;
    spheres[0].center[2] = 0.0;
    spheres[0].radius = 32.0;
    spheres[0].color[0] = 0.0;
    spheres[0].color[1] = 1.0;
    spheres[0].color[2] = 0.0;
}

void initTriangles() {
    triangles[0].p0[0] = 0;
    triangles[0].p0[1] = 25;
    triangles[0].p0[2] = 30;
    triangles[0].p1[0] = 32;
    triangles[0].p1[1] = 50;
    triangles[0].p1[2] = 20;
    triangles[0].p2[0] = -32;
    triangles[0].p2[1] = 50;
    triangles[0].p2[2] = 20;
    triangles[0].color[0] = 1.0;
    triangles[0].color[1] = 0.0;
    triangles[0].color[2] = 0.0;
    float p0[3] = {triangles[0].p0[0], triangles[0].p0[1], triangles[0].p0[2]};
    float p1[3] = {triangles[0].p1[0], triangles[0].p1[1], triangles[0].p1[2]};
    float p2[3] = {triangles[0].p2[0], triangles[0].p2[1], triangles[0].p2[2]};
    e1[0] = p1[0] - p0[0];
    e1[1] = p1[1] - p0[1];
    e1[2] = p1[2] - p0[2];
    
    e2[0] = p2[0] - p1[0];
    e2[1] = p2[1] - p1[1];
    e2[2] = p2[2] - p1[2];
    
    e3[0] = p2[0] - p0[0];
    e3[1] = p2[1] - p0[1];
    e3[2] = p2[2] - p0[2];
}

void drawPixels() {
    Ray ray;
    ray.p[0] = 0;
    ray.p[1] = -1000;
    ray.p[2] = 0;
    
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    for(float i = -128; i < 128; i++) {
        for(float j = -128; j < 128; j++) {
            ray.d[0] = i - ray.p[0];
            ray.d[1] = 0 - ray.p[1];
            ray.d[2] = j - ray.p[2];
            
            normalize(ray.d);
            
            rayTrace(ray.p, ray.d);
            glBegin(GL_POINTS);
            glColor3f(pix_col[0], pix_col[1], pix_col[2]);
            glVertex3f(i, j, 0);
            glEnd();
        }
    }
    
    glutSwapBuffers();
}

/** RAY TRACING FUNCTIONS **/

void rayTrace(float p[], float d[]) {
    float x[3];
    float n[3];
    
    float t_sph = findSphereIntersection(p, d);
    
    float t_tri = findTriangleIntersection(p, d);
    
    if(t_tri < 0.0) {
        if(t_sph > 0.0) {
            x[0] = p[0] + t_sph * d[0];
            x[1] = p[1] + t_sph * d[1];
            x[2] = p[2] + t_sph * d[2];
            
            n[0] = x[0] - spheres[0].center[0];
            n[1] = x[1] - spheres[1].center[1];
            n[2] = x[2] - spheres[2].center[2];
            
            normalize(n);
            
            pix_col[0] = spheres[0].color[0];
            pix_col[1] = spheres[0].color[1];
            pix_col[2] = spheres[0].color[2];
            
            calculateLighting(p, x, n, spheres[0].color);
        }
        else {
            pix_col[0] = 0.0;
            pix_col[1] = 0.0;
            pix_col[2] = 0.0;
        }
    }
    else {
        pix_col[0] = triangles[0].color[0];
        pix_col[1] = triangles[0].color[1];
        pix_col[2] = triangles[0].color[2];
        
        float n[3];
        crossProduct(e1, e2, n);
        normalize(n);
        x[0] = p[0] + t_tri * d[0];
        x[1] = p[1] + t_tri * d[1];
        x[2] = p[2] + t_tri * d[2];
        
        calculateLighting(p, x, n, triangles[0].color);
    }
}

float findSphereIntersection(float p[], float d[]) {
    float t1;
    float t2;
    float m[3];
    m[0] = p[0] - spheres[0].center[0];
    m[1] = p[1] - spheres[0].center[1];
    m[2] = p[2] - spheres[0].center[2];
    
    float r = spheres[0].radius;
    
    float discriminant = pow(dotProduct(d, m), 2) - (dotProduct(m, m) - pow(r, 2));
    t1 = -1 * (dotProduct(d, m)) + sqrt(discriminant);
    t2 = -1 * (dotProduct(d, m)) - sqrt(discriminant);
    if(discriminant > 0) {
        if(t1 < t2){
            return t1;
        }
        else {
            return t2;
        }
    }
    else if(discriminant == 0) {
        return t1;
    }
    else {
        return -1.0;
    }
}

float findTriangleIntersection(float p[], float d[]) {
    float p0[3] = {triangles[0].p0[0], triangles[0].p0[1], triangles[0].p0[2]};
    float p1[3] = {triangles[0].p1[0], triangles[0].p1[1], triangles[0].p1[2]};
    float p2[3] = {triangles[0].p2[0], triangles[0].p2[1], triangles[0].p2[2]};
    float normal[3];
    crossProduct(e1, e2, normal);

    normalize(normal);
    
    float k = dotProduct(normal, triangles[0].p1);
    
    float t = (k - dotProduct(p, normal)) / (dotProduct(d, normal));
    
    float r[3];
    
    r[0] = p[0] + t * d[0];
    r[1] = p[1] + t * d[1];
    r[2] = p[2] + t * d[2];
    
    float check1, check2, check3;
    
    float r1[3];
    r1[0] = r[0] - p0[0];
    r1[1] = r[1] - p0[1];
    r1[2] = r[2] - p0[2];
    
    float cross1[3];
    crossProduct(e1, r1, cross1);

    check1 = dotProduct(cross1, normal);
    
    /**/
    
    float r2[3];
    r2[0] = r[0] - p1[0];
    r2[1] = r[1] - p1[1];
    r2[2] = r[2] - p1[2];
    
    float cross2[3];
    crossProduct(e2, r2, cross2);
    
    check2 = dotProduct(cross2, normal);
    
    /**/
    
    float r3[3];
    r3[0] = r[0] - p2[0];
    r3[1] = r[1] - p2[1];
    r3[2] = r[2] - p2[2];
    
    float cross3[3];
    crossProduct(e3, r3, cross3);
    
    check3 = dotProduct(cross3, normal);
    
    if(check1 > 0 && check2 > 0 && check3 > 0) {
        return t;
    }
    else {
        return -1.0;
    }
}

void calculateLighting(float p[], float x[], float n[], float mat_col[]) {
    float l[3];
    l[0] = light_pos[0] - x[0];
    l[1] = light_pos[1] - x[1];
    l[2] = light_pos[2] - x[2];
    normalize(l);
    float amb_c[3];
    for(int i = 0; i < 3; i++) {
        amb_c[i] = light_color[i] * mat_col[i];
    }
    
    float diffuse_c[3];
    for(int i = 0; i < 3; i++) {
        diffuse_c[i] = light_color[i] * mat_col[i] * dotProduct(n, l);
    }
    
    float v[3];
    v[0] = p[0] - x[0];
    v[1] = p[1] - x[1];
    v[2] = p[2] - x[2];
    normalize(v);
    
    float h[3];
    h[0] = l[0] + v[0];
    h[1] = l[1] + v[1];
    h[2] = l[2] + v[2];
    normalize(h);
    
    float spec_c[3];
    for(int i = 0; i < 3; i++) {
        spec_c[i] = light_color[i] * mat_col[i] * dotProduct(n, h);
    }
    
    float final_color[3];
    final_color[0] = amb_c[0] + diffuse_c[0] + spec_c[0];
    final_color[1] = amb_c[1] + diffuse_c[1] + spec_c[1];
    final_color[2] = amb_c[2] + diffuse_c[2] + spec_c[2];
    
    pix_col[0] = final_color[0];
    pix_col[1] = final_color[1];
    pix_col[2] = final_color[2];
}

/** UTILITY FUNCTIONS **/

float dotProduct(float a[], float b[]) {
    float product = 0.0;
    for(int i = 0; i < 3; i++) {
        product += a[i] * b[i];
    }
    return product;
}

void crossProduct(float a[], float b[], float (&c)[3]) {
     c[0] = (a[1] * b[2]) - (a[2] * b[1]);
     c[1] = (a[2] * b[0]) - (a[0] * b[2]);
     c[2] = (a[0] * b[1]) - (a[1] * b[0]);
}

float findMagnitude(float a, float b, float c) {
    float mag = sqrt(pow(a, 2) + pow(b, 2) + pow(c, 2));
    return mag;
}

void normalize(float (&v)[3]) {
    float mag = findMagnitude(v[0], v[1], v[2]);
    v[0] /= mag;
    v[1] /= mag;
    v[2] /= mag;
}
