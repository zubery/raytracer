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

struct Vector3
{
  float x;
  float y;
  float z;
};

struct Color
{
  char r;
  char g;
  char b;
};

typedef struct _Triangle
{
  struct Vertex v[3];
} Triangle;

typedef struct _Sphere
{
  double position[3];
  double color_diffuse[3];
  double color_specular[3];
  double shininess;
  double radius;
} Sphere;

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

void create_frustum();
struct Vector3 lerpRay(int x, int y);
float dotProduct(struct Vector3 v1, struct Vector3 v2);
struct Vector3 crossProduct(struct Vector3 v1, struct Vector3 v2);
struct Vector3 segment(struct Vector3 point1, struct Vector3 point2);
struct Vector3 normalize(struct Vector3 original); 
void triangle_intersection(struct Vector3 direction, struct Hit& hit);
void sphere_intersection(struct Vector3 direction, struct Hit& hit);
bool pointInTriangle(struct Vector3 point, struct _Triangle tri); 
struct Color PhongSphere(struct Vector3 point, struct _Sphere sph); 

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
      struct Vector3 direction; 
      struct Color color;

      struct Hit triHit; 
      triHit.found = false; 
      triHit.t = std::numeric_limits<float>::max(); 
      struct Hit sphHit; 
      sphHit.found = false; 
      sphHit.t = std::numeric_limits<float>::max(); 

      direction = lerpRay(x, y);

      color.r = 255;
      color.g = 255;
      color.b = 255;

      triangle_intersection(direction, triHit);
      sphere_intersection(direction, sphHit); 

      if(triHit.found && triHit.t < sphHit.t)
      {
        color.r = 0; 
        color.g = 0; 
        color.b = 0; 
      }
      else if(sphHit.found && sphHit.t < triHit.t)
      {
        struct Vector3 point; 
        point.x = camera.x + (direction.x * sphHit.t); 
        point.y = camera.y + (direction.y * sphHit.t); 
        point.z = camera.z + (direction.z * sphHit.t); 

        color = PhongSphere(point, sphHit.sph); 
      }

      plot_pixel(x,y,color.r,color.g,color.b);
    }
    glEnd();
    glFlush();
  }
  printf("Done!\n"); fflush(stdout);
}

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

struct Vector3 lerpRay(int x, int y)
{
  struct Vector3 result;

  float xratio = ((float)x) / ((float)WIDTH);
  float yratio = ((float)y) / ((float)HEIGHT);

  struct Vector3 left;
  left.x = (ray00.x + ray01.x) / 2.0f; 
  left.y = ray00.y + (yratio * (ray01.y - ray00.y));
  left.z = -1.0f;

  struct Vector3 right;
  right.x = (ray10.x + ray11.x) / 2.0f;
  right.y = ray10.y + (yratio * (ray11.y - ray10.y));
  right.z = -1.0f;

  result.x = left.x + (xratio * (right.x - left.x));
  result.y = (left.y + right.y) / 2.0f;
  result.z = -1.0f;

  float magnitude = sqrt(pow(result.x, 2.0f) + pow(result.y, 2.0f) + pow(result.z, 2.0f));

  result.x = result.x / magnitude; 
  result.y = result.y / magnitude; 
  result.z = result.z / magnitude; 

  return result; 
}

float dotProduct(struct Vector3 v1, struct Vector3 v2)
{
  float result = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z);
  return result; 
}

struct Vector3 crossProduct(struct Vector3 v1, struct Vector3 v2)
{
  struct Vector3 result; 

  result.x = (v1.y * v2.z) - (v1.z * v2.y); 
  result.y = (v1.z * v2.x) - (v1.x * v2.z);
  result.z = (v1.x * v2.y) - (v1.y * v2.x); 

  return result; 
}

struct Vector3 segment(struct Vector3 point1, struct Vector3 point2)
{
  struct Vector3 result; 

  result.x = point2.x - point1.x; 
  result.y = point2.y - point1.y; 
  result.z = point2.z - point1.z; 

  return result; 
}

struct Vector3 normalize(struct Vector3 original)
{
  float length = sqrt(dotProduct(original, original)); 
  
  struct Vector3 result; 
  result.x = original.x / length; 
  result.y = original.y / length; 
  result.z = original.z / length; 

  return result; 
}

void triangle_intersection(struct Vector3 direction, struct Hit& hit)
{
  for(int i = 0; i < num_triangles; i++)
  {
    struct _Triangle current = triangles[i];

    struct Vector3 normal0; 
    normal0.x = current.v[0].normal[0]; 
    normal0.y = current.v[0].normal[1]; 
    normal0.z = current.v[0].normal[2]; 
    
    struct Vector3 vertex0;
    vertex0.x = current.v[0].position[0];
    vertex0.y = current.v[0].position[1];
    vertex0.z = current.v[0].position[2];

    float t; 
    float top; 
    float bottom = dotProduct(normal0, direction);

    if(!(fabs(0.0f - bottom) < 0.005f))
    {
      float d = dotProduct(normal0, vertex0); 
      top = (dotProduct(normal0, camera) + d); 

      t = top / bottom; 

      if(t >= 0)
      {
        struct Vector3 intersection; 
        intersection.x = camera.x + (direction.x * t); 
        intersection.y = camera.y + (direction.y * t); 
        intersection.z = camera.z + (direction.z * t); 

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

void sphere_intersection(struct Vector3 direction, struct Hit& hit)
{
  for(int i = 0; i < num_spheres; i++)
  {
    struct _Sphere current = spheres[i];

    struct Vector3 center; 
    center.x = current.position[0]; 
    center.y = current.position[1];
    center.z = current.position[2];

    float a = pow(direction.x, 2.0f) + pow(direction.y, 2.0f) + pow(direction.z, 2.0f); 
    float b = 2.0f * ((direction.x * (camera.x - center.x)) + (direction.y * (camera.y - center.y)) + (direction.z * (camera.z - center.z))); 
    float c = pow(camera.x - center.x, 2.0f) + pow(camera.y - center.y, 2.0f) + pow(camera.z - center.z, 2.0f) - pow(current.radius, 2.0f);

    float delta = pow(b, 2.0f) - (4 * a * c); 
    float t0; 
    float t1;

    if(delta >= 0)
    {
      t0 = ((-1 * b) + (sqrt(delta))) / 2.0f;
      t1 = ((-1 * b) - (sqrt(delta))) / 2.0f;

      if(hit.found == false)
      {
        hit.found = true; 
        hit.t = std::min(t0, t1);
        hit.sph = current; 
      }
      else
      {
        float tmin = std::min(t0, t1); 
        if(tmin < hit.t)
        {
          hit.t = tmin; 
          hit.sph = current; 
        }
      }
    }
  }
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
  else if(phongG > 1.0f)
  {
    phongG = 1.0f; 
  }
  else if(phongB > 1.0f)
  {
    phongB = 1.0f; 
  }

  struct Color result;

  result.r = (int)(255 * phongR); 
  result.g = (int)(255 * phongR); 
  result.b = (int)(255 * phongR);

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
