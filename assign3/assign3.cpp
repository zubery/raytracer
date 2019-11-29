/*
CSCI 480
Assignment 3 Raytracer

Name: Ryan Zubery
*/

#include <stdlib.h>
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>
#include <GLUT/glut.h>
#include <pic.h>
#include <string.h>
#include <math.h>
#include <iostream>
#include <algorithm>
#include <limits>

#define MAX_TRIANGLES 2000
#define MAX_SPHERES 10
#define MAX_LIGHTS 10

char *filename=0;

//different display modes
#define MODE_DISPLAY 1
#define MODE_JPEG 2
int mode=MODE_DISPLAY;

//you may want to make these smaller for debugging purposes
#define WIDTH 640
#define HEIGHT 480

//the field of view of the camera
#define fov 60.0

unsigned char buffer[HEIGHT][WIDTH][3];

struct Vertex
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double normal[3];
  double shininess;
};


//3 points
struct Vector3
{
  float x;
  float y;
  float z;
};

//container for rgb color
struct Color
{
  char r;
  char g;
  char b;
};

typedef struct _Triangle
{
  struct Vertex v[3];
  int id; 
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
  int id; 
} Sphere;

//container for raycast hit information
struct Hit
{
  bool found; 
  float t; 
  struct _Triangle tri; 
  struct _Sphere sph; 
}; 

typedef struct _Light
{
  double position[3];
  double color[3];
} Light;

Triangle triangles[MAX_TRIANGLES];
Sphere spheres[MAX_SPHERES];
Light lights[MAX_LIGHTS];
double ambient_light[3];

int num_triangles=0;
int num_spheres=0;
int num_lights=0;

struct Vector3 camera;
struct Vector3 ray00;
struct Vector3 ray01;
struct Vector3 ray10;
struct Vector3 ray11;

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b);
void plot_pixel(int x,int y,unsigned char r,unsigned char g,unsigned char b);

//functions to create raycasts
void create_frustum();
struct Vector3 lerpRay(int x, int y);

//linear algebra and 3d mathfunctions
float dotProduct(struct Vector3 v1, struct Vector3 v2);
struct Vector3 crossProduct(struct Vector3 v1, struct Vector3 v2);
struct Vector3 segment(struct Vector3 point1, struct Vector3 point2);
struct Vector3 normalize(struct Vector3 original); 

//intersection functions
void triangle_intersection(struct Vector3 direction, struct Vector3 origin, struct Hit& hit, int ignoreID);
void sphere_intersection(struct Vector3 direction, struct Vector3 origin, struct Hit& hit, int ignoreID);
bool pointInTriangle(struct Vector3 point, struct _Triangle tri); 

//barycentric functions
struct Vector3 baryWeights(struct Vector3 point, struct Vector3 v1, struct Vector3 v2, struct Vector3 v3); 
struct Vector3 interpolateTriangle(struct Vector3 val1, struct Vector3 val2, struct Vector3 val3, struct Vector3 weights); 

//phong shading functions
struct Color PhongSphere(struct Vector3 point, struct _Sphere sph); 
struct Color PhongTriangle(struct Vector3 point, struct _Triangle tri);

//MODIFY THIS FUNCTION
void draw_scene()
{
  unsigned int x,y;

  create_frustum();

  //simple output
  for(x=0; x<WIDTH; x++)
  {
    glPointSize(2.0);  
    glBegin(GL_POINTS);
    for(y=0;y < HEIGHT;y++)
    {
      //initializing variables for raycasting
      struct Vector3 direction; 
      struct Color color;

      struct Hit triHit; 
      triHit.found = false; 
      triHit.t = std::numeric_limits<float>::max(); 
      struct Hit sphHit; 
      sphHit.found = false; 
      sphHit.t = std::numeric_limits<float>::max(); 

      //creates the raycast for this pixel
      direction = lerpRay(x, y);

      //default color value if no raycast hit
      color.r = 255;
      color.g = 255;
      color.b = 255;

      //triangle and sphere intersection tests
      triangle_intersection(direction, camera, triHit, -1);
      sphere_intersection(direction, camera, sphHit, -1); 

      //triangle found, intersection point shaded as triangle
      if(triHit.found && triHit.t < sphHit.t)
      {
        struct Vector3 point; 
        point.x = camera.x + (direction.x * triHit.t); 
        point.y = camera.y + (direction.y * triHit.t); 
        point.z = camera.z + (direction.z * triHit.t); 

        color = PhongTriangle(point, triHit.tri); 
      }
      //triangle found, intersection point shaded as sphere
      else if(sphHit.found && sphHit.t < triHit.t)
      {
        struct Vector3 point; 
        point.x = camera.x + (direction.x * sphHit.t); 
        point.y = camera.y + (direction.y * sphHit.t); 
        point.z = camera.z + (direction.z * sphHit.t); 

        color = PhongSphere(point, sphHit.sph); 
      }

      //pixel is plotted
      plot_pixel(x,y,color.r,color.g,color.b);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

//creates the 4 rays that bound the viewing window
void create_frustum()
{
  camera.x = 0.0f;
  camera.y = 0.0f;
  camera.z = 0.0f;

  float aspect = ((float)WIDTH) / ((float)HEIGHT);

  float halfFOV = fov / 2.0f;
  float halfFOVRad = halfFOV * (M_PI / 180.0f);
  float tanHalfFOVRad = tan(halfFOVRad);

  ray00.x = -1.0f * aspect * tanHalfFOVRad;
  ray00.y = -1.0f * tanHalfFOVRad;
  ray00.z = -1.0f;

  ray01.x = -1.0f * aspect * tanHalfFOVRad;
  ray01.y = tanHalfFOVRad;
  ray01.z = -1.0f;

  ray10.x = aspect * tanHalfFOVRad;
  ray10.y = -1.0f * tanHalfFOVRad;
  ray10.z = -1.0f;

  ray11.x = aspect * tanHalfFOVRad;
  ray11.y = tanHalfFOVRad;
  ray11.z = -1.0f;
}

//for a given pixel in the window, lerps between the 4 bounding rays for specific pixel ray
struct Vector3 lerpRay(int x, int y)
{
  struct Vector3 result;

  float xratio = ((float)x) / ((float)WIDTH);
  float yratio = ((float)y) / ((float)HEIGHT);

  //lerp first in terms of height
  struct Vector3 left;
  left.x = (ray00.x + ray01.x) / 2.0f; 
  left.y = ray00.y + (yratio * (ray01.y - ray00.y));
  left.z = -1.0f;

  struct Vector3 right;
  right.x = (ray10.x + ray11.x) / 2.0f;
  right.y = ray10.y + (yratio * (ray11.y - ray10.y));
  right.z = -1.0f;

  //lerp in terms of width
  result.x = left.x + (xratio * (right.x - left.x));
  result.y = (left.y + right.y) / 2.0f;
  result.z = -1.0f;

  float magnitude = sqrt(pow(result.x, 2.0f) + pow(result.y, 2.0f) + pow(result.z, 2.0f));

  result.x = result.x / magnitude; 
  result.y = result.y / magnitude; 
  result.z = result.z / magnitude; 

  return result; 
}

//dot product of two 3d vectors
float dotProduct(struct Vector3 v1, struct Vector3 v2)
{
  float result = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
  return result; 
}

//cross product of two 3d vectors
struct Vector3 crossProduct(struct Vector3 v1, struct Vector3 v2)
{
  struct Vector3 result; 

  result.x = (v1.y * v2.z) - (v1.z * v2.y); 
  result.y = (v1.z * v2.x) - (v1.x * v2.z);
  result.z = (v1.x * v2.y) - (v1.y * v2.x); 

  return result; 
}

//vector between two 3d points
struct Vector3 segment(struct Vector3 point1, struct Vector3 point2)
{
  struct Vector3 result; 

  result.x = point2.x - point1.x; 
  result.y = point2.y - point1.y; 
  result.z = point2.z - point1.z; 

  return result; 
}

//normalized vector
struct Vector3 normalize(struct Vector3 original)
{
  float length = sqrt(dotProduct(original, original)); 
  
  struct Vector3 result; 
  result.x = original.x / length; 
  result.y = original.y / length; 
  result.z = original.z / length; 

  return result; 
}

//checks for ray intersection with triangle
void triangle_intersection(struct Vector3 direction, struct Vector3 origin, struct Hit& hit, int ignoreID)
{
  //checks each triangle in the scene
  for(int i = 0; i < num_triangles; i++)
  {
    struct _Triangle current = triangles[i];
    
    if(current.id == ignoreID)
    {
      continue; 
    }

    struct Vector3 vertex0;
    vertex0.x = current.v[0].position[0];
    vertex0.y = current.v[0].position[1];
    vertex0.z = current.v[0].position[2];

    struct Vector3 vertex1;
    vertex1.x = current.v[1].position[0];
    vertex1.y = current.v[1].position[1];
    vertex1.z = current.v[1].position[2];

    struct Vector3 vertex2;
    vertex2.x = current.v[2].position[0];
    vertex2.y = current.v[2].position[1];
    vertex2.z = current.v[2].position[2];

    struct Vector3 edge01 = segment(vertex0, vertex1); 
    struct Vector3 edge02 = segment(vertex0, vertex2); 

    //calculates normal
    struct Vector3 normal0 = crossProduct(edge01, edge02); 
    normal0 = normalize(normal0); 

    //solves for intersection point
    float t; 
    float top; 
    float bottom = dotProduct(normal0, direction);

    if(!(fabs(0.0f - bottom) < 0.005f))
    {
      float d = -1.0f * dotProduct(normal0, vertex0); 
      top = -1.0f * (dotProduct(normal0, origin) + d); 

      t = top / bottom; 

      if(t >= 0)
      {
        struct Vector3 intersection; 
        intersection.x = origin.x + (direction.x * t); 
        intersection.y = origin.y + (direction.y * t); 
        intersection.z = origin.z + (direction.z * t); 

        //std::cout << intersection.x << ", " << intersection.y << ", " << intersection.z << std::endl; 

        bool inside = pointInTriangle(intersection, current); 

        if(inside)
        {
          if(hit.found == false)
          {
            hit.found = true; 
            hit.t = t;
            hit.tri = current; 
          }
          else if(t < hit.t)
          {
            hit.t = t; 
            hit.tri = current; 
          }
        }
      }
    }
  }
}

//checks if point in triangle by using barycentric coordinates and linear algebra
bool pointInTriangle(struct Vector3 point, struct _Triangle tri)
{
  struct Vector3 vertA; 
  vertA.x = tri.v[0].position[0]; 
  vertA.y = tri.v[0].position[1]; 
  vertA.z = tri.v[0].position[2]; 

  struct Vector3 vertB; 
  vertB.x = tri.v[1].position[0]; 
  vertB.y = tri.v[1].position[1]; 
  vertB.z = tri.v[1].position[2]; 

  struct Vector3 vertC; 
  vertC.x = tri.v[2].position[0]; 
  vertC.y = tri.v[2].position[1]; 
  vertC.z = tri.v[2].position[2]; 

  struct Vector3 u = segment(vertA, vertB); 
  struct Vector3 v = segment(vertA, vertC); 
  struct Vector3 w = segment(vertA, point); 

  struct Vector3 vCrossW = crossProduct(v, w); 
  struct Vector3 vCrossU = crossProduct(v, u); 

  if(dotProduct(vCrossW, vCrossU) < 0.0f)
  {
    return false; 
  }

  struct Vector3 uCrossW = crossProduct(u, w); 
  struct Vector3 uCrossV = crossProduct(u, v); 

  if(dotProduct(uCrossW, uCrossV) < 0.0f)
  {
    return false; 
  }

  float divideBy = sqrt(dotProduct(uCrossV, uCrossV)); 
  float r = sqrt(dotProduct(vCrossW, vCrossW)) / divideBy; 
  float t = sqrt(dotProduct(uCrossW, uCrossW)) / divideBy; 

  return (r + t <= 1); 
}

//returns respective weights of barycentric coordinates
struct Vector3 baryWeights(struct Vector3 point, struct Vector3 v1, struct Vector3 v2, struct Vector3 v3)
{
  struct Vector3 weight; 

  weight.x = ((v2.y - v3.y)*(point.x - v3.x) + (v3.x - v2.x)*(point.y - v3.y)) / ((v2.y - v3.y)*(v1.x - v3.x) + (v3.x - v2.x)*(v1.y - v3.y)); 

  weight.y = ((v3.y - v1.y)*(point.x - v3.x) + (v1.x - v3.x)*(point.y - v3.y)) / ((v2.y - v3.y)*(v1.x - v3.x) + (v3.x - v2.x)*(v1.y - v3.y)); 

  weight.z = 1 - weight.x - weight.y; 

  return weight; 
}

//interpolates a point on a triangle given the vertices and weights
struct Vector3 interpolateTriangle(struct Vector3 val1, struct Vector3 val2, struct Vector3 val3, struct Vector3 weights)
{
  float weight1 = weights.x; 
  float weight2 = weights.y; 
  float weight3 = weights.z; 

  if(weight1 < 0.0f)
  {
    weight1 = 0.0f; 
  }
  else if(weight2 < 0.0f)
  {
    weight2 = 0.0f; 
  }
  else if(weight3 < 0.0f)
  {
    weight3 = 0.0f; 
  }

  struct Vector3 results; 
  results.x = (val1.x * weight1) + (val2.x * weight2) + (val3.x * weight3); 
  results.y = (val1.y * weight1) + (val2.y * weight2) + (val3.y * weight3); 
  results.z = (val1.z * weight1) + (val2.z * weight2) + (val3.z * weight3); 

  return results; 
}

//calculates sphere intersection
void sphere_intersection(struct Vector3 direction, struct Vector3 origin, struct Hit& hit, int ignoreID)
{
  //checks each sphere in the scene
  for(int i = 0; i < num_spheres; i++)
  {
    struct _Sphere current = spheres[i];

    if(current.id == ignoreID)
    {
      continue; 
    }

    struct Vector3 center; 
    center.x = current.position[0]; 
    center.y = current.position[1];
    center.z = current.position[2];

    float a = pow(direction.x, 2.0f) + pow(direction.y, 2.0f) + pow(direction.z, 2.0f); 
    float b = 2.0f * ((direction.x * (origin.x - center.x)) + (direction.y * (origin.y - center.y)) + (direction.z * (origin.z - center.z))); 
    float c = pow(origin.x - center.x, 2.0f) + pow(origin.y - center.y, 2.0f) + pow(origin.z - center.z, 2.0f) - pow(current.radius, 2.0f);

    float delta = pow(b, 2.0f) - (4 * a * c); 
    float t0; 
    float t1;

    //calculates intersection point
    if(delta >= 0)
    {
      t0 = ((-1 * b) + (sqrt(delta))) / 2.0f;
      t1 = ((-1 * b) - (sqrt(delta))) / 2.0f;

      bool check = (t0 < 0) && (t1 < 0); 

      if(!check)
      {
        if(hit.found == false)
        {
          hit.t = std::min(t0, t1);
          
          if(hit.t < 0)
          {
            hit.t = std::max(t0, t1); 
          }

          hit.found = true; 
          hit.sph = current; 
        }
        else
        {
          float tmin = std::min(t0, t1); 

          if(tmin < 0)
          {
            tmin = std::max(t0, t1);
          }

          if(tmin < hit.t)
          {
            hit.t = tmin; 
            hit.sph = current; 
          }
        }
      }
    }
  }
}

//phong shading for points in triangle
struct Color PhongTriangle(struct Vector3 point, struct _Triangle tri)
{
  //VERTICES POSITIONS
  struct Vector3 v1; 
  v1.x = tri.v[0].position[0]; 
  v1.y = tri.v[0].position[1]; 
  v1.z = tri.v[0].position[2]; 

  struct Vector3 v2; 
  v2.x = tri.v[1].position[0]; 
  v2.y = tri.v[1].position[1]; 
  v2.z = tri.v[1].position[2]; 

  struct Vector3 v3; 
  v3.x = tri.v[2].position[0]; 
  v3.y = tri.v[2].position[1]; 
  v3.z = tri.v[2].position[2]; 

  //FIND THE WEIGHTS
  struct Vector3 weights = baryWeights(point, v1, v2, v3); 

  //VERT 1 PROPERTIES
  struct Vector3 v1Diffuse; 
  v1Diffuse.x = tri.v[0].color_diffuse[0]; 
  v1Diffuse.y = tri.v[0].color_diffuse[1]; 
  v1Diffuse.z = tri.v[0].color_diffuse[2]; 

  struct Vector3 v1Specular;
  v1Specular.x = tri.v[0].color_specular[0]; 
  v1Specular.y = tri.v[0].color_specular[1]; 
  v1Specular.z = tri.v[0].color_specular[2]; 

  struct Vector3 v1Normal; 
  v1Normal.x = tri.v[0].normal[0]; 
  v1Normal.y = tri.v[0].normal[1]; 
  v1Normal.z = tri.v[0].normal[2];

  float v1Shine = tri.v[0].shininess;  

  //VERT 2 PROPERTIES
  struct Vector3 v2Diffuse; 
  v2Diffuse.x = tri.v[1].color_diffuse[0]; 
  v2Diffuse.y = tri.v[1].color_diffuse[1]; 
  v2Diffuse.z = tri.v[1].color_diffuse[2]; 

  struct Vector3 v2Specular;
  v2Specular.x = tri.v[1].color_specular[0]; 
  v2Specular.y = tri.v[1].color_specular[1]; 
  v2Specular.z = tri.v[1].color_specular[2]; 

  struct Vector3 v2Normal; 
  v2Normal.x = tri.v[1].normal[0]; 
  v2Normal.y = tri.v[1].normal[1]; 
  v2Normal.z = tri.v[1].normal[2];

  float v2Shine = tri.v[1].shininess; 

  //VERT 3 PROPERTIES
  struct Vector3 v3Diffuse; 
  v3Diffuse.x = tri.v[2].color_diffuse[0]; 
  v3Diffuse.y = tri.v[2].color_diffuse[1]; 
  v3Diffuse.z = tri.v[2].color_diffuse[2]; 

  struct Vector3 v3Specular;
  v3Specular.x = tri.v[2].color_specular[0]; 
  v3Specular.y = tri.v[2].color_specular[1]; 
  v3Specular.z = tri.v[2].color_specular[2]; 

  struct Vector3 v3Normal; 
  v3Normal.x = tri.v[2].normal[0]; 
  v3Normal.y = tri.v[2].normal[1]; 
  v3Normal.z = tri.v[2].normal[2];

  float v3Shine = tri.v[2].shininess; 

  //INTERPOLATED VALUES
  struct Vector3 interDiff = interpolateTriangle(v1Diffuse, v2Diffuse, v3Diffuse, weights); 

  struct Vector3 interSpec = interpolateTriangle(v1Specular, v2Specular, v3Specular, weights); 

  struct Vector3 interNorm = interpolateTriangle(v1Normal, v2Normal, v3Normal, weights); 

  //RGB VALUES
  float kSpecularR = interSpec.x; 
  float kSpecularG = interSpec.y; 
  float kSpecularB = interSpec.z; 

  float kDiffuseR = interDiff.x; 
  float kDiffuseG = interDiff.y; 
  float kDiffuseB = interDiff.z; 

  float kAmbientR = ambient_light[0]; 
  float kAmbientG = ambient_light[1]; 
  float kAmbientB = ambient_light[2]; 

  //SHINE CALCULATION
  float w1 = weights.x; 
  float w2 = weights.y; 
  float w3 = weights.z; 

  if(w1 < 0.0f)
  {
    w1 = 0.0f; 
  }
  else if(w2 < 0.0f)
  {
    w2 = 0.0f; 
  }
  else if(w3 < 0.0f)
  {
    w3 = 0.0f; 
  }

  float shine = (v1Shine * w1) + (v2Shine * w2) + (v3Shine * w3); 

  //PHONG CALCULATION, SIMILAR TO BELOW SPHERE CALCULATIONS
  float phongR = kAmbientR; 
  float phongG = kAmbientG; 
  float phongB = kAmbientB; 

  for(int i = 0; i < num_lights; i++)
  {
    struct _Light current = lights[i];

    struct Vector3 lPoint; 
    lPoint.x = current.position[0]; 
    lPoint.y = current.position[1]; 
    lPoint.z = current.position[2]; 

    float currR = current.color[0]; 
    float currG = current.color[1]; 
    float currB = current.color[2]; 

    //CASTS TO SEE IF THERES A SHADOW BETWEEN THIS LIGHT
    struct Vector3 shadow = segment(point, lPoint); 
    shadow = normalize(shadow); 
    struct Hit triHit; 
    triHit.found = false; 
    struct Hit sphHit; 
    sphHit.found = false; 
    triangle_intersection(shadow, point, triHit, tri.id); 
    sphere_intersection(shadow, point, sphHit, -1); 

    if(sphHit.found)
    {
      currR = 0.0f; 
      currG = 0.0f; 
      currB = 0.0f; 
    }

    if(triHit.found)
    {
      if(triHit.tri.id != tri.id)
      {
        currR = 0.0f; 
        currG = 0.0f; 
        currB = 0.0f;
      }
    }
    //DONE WITH SHADOW CAST

    struct Vector3 L = segment(point, lPoint); 
    L = normalize(L);  

    struct Vector3 normal = interNorm; 
    normal = normalize(normal); 

    struct Vector3 REnd; 
    REnd.x = 2 * dotProduct(normal, L) * normal.x; 
    REnd.y = 2 * dotProduct(normal, L) * normal.y; 
    REnd.z = 2 * dotProduct(normal, L) * normal.z; 

    struct Vector3 R = segment(L, REnd); 
    R = normalize(R); 

    struct Vector3 V = segment(point, camera); 
    V = normalize(V); 

    float lDotN = dotProduct(L, normal); 
    if(lDotN < 0.0f)
    {
      lDotN = 0.0f; 
    }

    float rDotV = dotProduct(R, V);
    if(rDotV < 0.0f)
    {
      rDotV = 0.0f; 
    }

    currR *= ((kDiffuseR * lDotN) + (kSpecularR * pow(rDotV, shine))); 

    currG *= ((kDiffuseG * lDotN) + (kSpecularG * pow(rDotV, shine))); 

    currB *= ((kDiffuseB * lDotN) + (kSpecularB * pow(rDotV, shine))); 

    phongR += currR; 
    phongG += currG; 
    phongB += currB; 
  }

  if(phongR > 1.0f)
  {
    phongR = 1.0f; 
  }
  if(phongG > 1.0f)
  {
    phongG = 1.0f; 
  }
  if(phongB > 1.0f)
  {
    phongB = 1.0f; 
  }

  struct Color result;

  result.r = (int)(255 * phongR); 
  result.g = (int)(255 * phongG); 
  result.b = (int)(255 * phongB);

  return result; 
}

struct Color PhongSphere(struct Vector3 point, struct _Sphere sph)
{
  float kSpecularR = sph.color_specular[0];
  float kDiffuseR = sph.color_diffuse[0]; 
  float kAmbientR = ambient_light[0]; 

  float kSpecularG = sph.color_specular[1];
  float kDiffuseG = sph.color_diffuse[1]; 
  float kAmbientG = ambient_light[1];

  float kSpecularB = sph.color_specular[2];
  float kDiffuseB = sph.color_diffuse[2]; 
  float kAmbientB = ambient_light[2];

  float shine = sph.shininess; 

  float phongR = kAmbientR; 
  float phongG = kAmbientG; 
  float phongB = kAmbientB; 

  struct Vector3 sphCenter; 
  sphCenter.x = sph.position[0]; 
  sphCenter.y = sph.position[1]; 
  sphCenter.z = sph.position[2];  

  //for each light, calculates contribution
  for(int i = 0; i < num_lights; i++)
  {
    struct _Light current = lights[i];

    struct Vector3 lPoint; 
    lPoint.x = current.position[0]; 
    lPoint.y = current.position[1]; 
    lPoint.z = current.position[2]; 

    float currR = current.color[0]; 
    float currG = current.color[1]; 
    float currB = current.color[2]; 

    //CASTS TO SEE IF THERES A SHADOW BETWEEN THIS LIGHT
    struct Vector3 shadow = segment(point, lPoint); 
    shadow = normalize(shadow); 

    struct Hit triHit; 
    triHit.found = false; 
    struct Hit sphHit; 
    sphHit.found = false; 
    triangle_intersection(shadow, point, triHit, -1); 
    sphere_intersection(shadow, point, sphHit, sph.id); 

    if(triHit.found)
    { 
      currR = 0.0f; 
      currG = 0.0f; 
      currB = 0.0f; 
    }

    if(sphHit.found)
    {
      if(sphHit.sph.id != sph.id)
      {
        currR = 0.0f; 
        currG = 0.0f; 
        currB = 0.0f;
      }
    }
    //DONE WITH SHADOW CAST

    struct Vector3 L = segment(point, lPoint); 
    L = normalize(L);  

    struct Vector3 normal = segment(sphCenter, point); 
    normal = normalize(normal); 

    struct Vector3 REnd; 
    REnd.x = 2 * dotProduct(normal, L) * normal.x; 
    REnd.y = 2 * dotProduct(normal, L) * normal.y; 
    REnd.z = 2 * dotProduct(normal, L) * normal.z; 

    struct Vector3 R = segment(L, REnd); 
    R = normalize(R); 

    struct Vector3 V = segment(point, camera); 
    V = normalize(V); 

    float lDotN = dotProduct(L, normal); 
    if(lDotN < 0.0f)
    {
      lDotN = 0.0f; 
    }

    float rDotV = dotProduct(R, V);
    if(rDotV < 0.0f)
    {
      rDotV = 0.0f; 
    }

    currR *= ((kDiffuseR * lDotN) + (kSpecularR * pow(rDotV, shine))); 

    currG *= ((kDiffuseG * lDotN) + (kSpecularG * pow(rDotV, shine))); 

    currB *= ((kDiffuseB * lDotN) + (kSpecularB * pow(rDotV, shine))); 

    phongR += currR; 
    phongG += currG; 
    phongB += currB; 
  }

  if(phongR > 1.0f)
  {
    phongR = 1.0f; 
  }
  if(phongG > 1.0f)
  {
    phongG = 1.0f; 
  }
  if(phongB > 1.0f)
  {
    phongB = 1.0f; 
  }

  struct Color result;

  result.r = (int)(255 * phongR); 
  result.g = (int)(255 * phongG); 
  result.b = (int)(255 * phongB);

  return result; 
}

void plot_pixel_display(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  glColor3f(((double)r)/256.f,((double)g)/256.f,((double)b)/256.f);
  glVertex2i(x,y);
}

void plot_pixel_jpeg(int x,int y,unsigned char r,unsigned char g,unsigned char b)
{
  buffer[HEIGHT-y-1][x][0]=r;
  buffer[HEIGHT-y-1][x][1]=g;
  buffer[HEIGHT-y-1][x][2]=b;
}

void plot_pixel(int x,int y,unsigned char r,unsigned char g, unsigned char b)
{
  plot_pixel_display(x,y,r,g,b);
  if(mode == MODE_JPEG)
      plot_pixel_jpeg(x,y,r,g,b);
}

void save_jpg()
{
  Pic *in = NULL;

  in = pic_alloc(640, 480, 3, NULL);
  printf("Saving JPEG file: %s\n", filename);

  memcpy(in->pix,buffer,3*WIDTH*HEIGHT);
  if (jpeg_write(filename, in))
    printf("File saved Successfully\n");
  else
    printf("Error in Saving\n");

  pic_free(in);      

}

void parse_check(char *expected,char *found)
{
  if(strcasecmp(expected,found))
    {
      char error[100];
      printf("Expected '%s ' found '%s '\n",expected,found);
      printf("Parse error, abnormal abortion\n");
      exit(0);
    }

}

void parse_doubles(FILE*file, char *check, double p[3])
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check(check,str);
  fscanf(file,"%lf %lf %lf",&p[0],&p[1],&p[2]);
  printf("%s %lf %lf %lf\n",check,p[0],p[1],p[2]);
}

void parse_rad(FILE*file,double *r)
{
  char str[100];
  fscanf(file,"%s",str);
  parse_check("rad:",str);
  fscanf(file,"%lf",r);
  printf("rad: %f\n",*r);
}

void parse_shi(FILE*file,double *shi)
{
  char s[100];
  fscanf(file,"%s",s);
  parse_check("shi:",s);
  fscanf(file,"%lf",shi);
  printf("shi: %f\n",*shi);
}

int loadScene(char *argv)
{
  FILE *file = fopen(argv,"r");
  int number_of_objects;
  char type[50];
  int i;
  Triangle t;
  Sphere s;
  Light l;
  fscanf(file,"%i",&number_of_objects);

  printf("number of objects: %i\n",number_of_objects);
  char str[200];

  parse_doubles(file,"amb:",ambient_light);

  for(i=0;i < number_of_objects;i++)
    {
      fscanf(file,"%s\n",type);
      printf("%s\n",type);
      if(strcasecmp(type,"triangle")==0)
	{

	  printf("found triangle\n");
	  int j;

	  for(j=0;j < 3;j++)
	    {
	      parse_doubles(file,"pos:",t.v[j].position);
	      parse_doubles(file,"nor:",t.v[j].normal);
	      parse_doubles(file,"dif:",t.v[j].color_diffuse);
	      parse_doubles(file,"spe:",t.v[j].color_specular);
	      parse_shi(file,&t.v[j].shininess);
	    }

	  if(num_triangles == MAX_TRIANGLES)
	    {
	      printf("too many triangles, you should increase MAX_TRIANGLES!\n");
	      exit(0);
	    }
	  t.id = num_triangles; 
    triangles[num_triangles++] = t;
	}
      else if(strcasecmp(type,"sphere")==0)
	{
	  printf("found sphere\n");

	  parse_doubles(file,"pos:",s.position);
	  parse_rad(file,&s.radius);
	  parse_doubles(file,"dif:",s.color_diffuse);
	  parse_doubles(file,"spe:",s.color_specular);
	  parse_shi(file,&s.shininess);

	  if(num_spheres == MAX_SPHERES)
	    {
	      printf("too many spheres, you should increase MAX_SPHERES!\n");
	      exit(0);
	    }
    s.id = num_spheres; 
	  spheres[num_spheres++] = s;
	}
      else if(strcasecmp(type,"light")==0)
	{
	  printf("found light\n");
	  parse_doubles(file,"pos:",l.position);
	  parse_doubles(file,"col:",l.color);

	  if(num_lights == MAX_LIGHTS)
	    {
	      printf("too many lights, you should increase MAX_LIGHTS!\n");
	      exit(0);
	    }
	  lights[num_lights++] = l;
	}
      else
	{
	  printf("unknown type in scene description:\n%s\n",type);
	  exit(0);
	}
    }
  return 0;
}

void display()
{

}

void reshape(int w, int h)
{
  //set up image size
  glViewport(0, 0, WIDTH * 2, HEIGHT * 2);
}

void init()
{
  glMatrixMode(GL_PROJECTION);
  glOrtho(0,WIDTH,0,HEIGHT,1,-1);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  glClearColor(0,0,0,0);
  glClear(GL_COLOR_BUFFER_BIT);
}

void idle()
{
  //hack to make it only draw once
  static int once=0;
  if(!once)
  {
      draw_scene();
      if(mode == MODE_JPEG)
	save_jpg();
    }
  once=1;
}

int main (int argc, char ** argv)
{
  if (argc<2 || argc > 3)
  {  
    printf ("usage: %s <scenefile> [jpegname]\n", argv[0]);
    exit(0);
  }
  if(argc == 3)
    {
      mode = MODE_JPEG;
      filename = argv[2];
    }
  else if(argc == 2)
    mode = MODE_DISPLAY;

  glutInit(&argc,argv);
  loadScene(argv[1]);

  glutInitDisplayMode(GLUT_RGBA | GLUT_SINGLE);
  glutInitWindowPosition(0,0);
  glutInitWindowSize(WIDTH,HEIGHT);
  int window = glutCreateWindow("Ray Tracer");
  glutDisplayFunc(display);
  glutReshapeFunc(reshape); 
  glutIdleFunc(idle);
  init();
  glutMainLoop();
}
