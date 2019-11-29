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
void triangle_intersection(struct Vector3 direction, struct Color color);
void sphere_intersection(struct Vector3 direction, struct Color& color);

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

      direction = lerpRay(x, y);

      color.r = 255;
      color.g = 255;
      color.b = 255;

      triangle_intersection(direction, color);
      sphere_intersection(direction, color); 

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

  float aspect = WIDTH/HEIGHT;
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

void triangle_intersection(struct Vector3 direction, struct Color color)
{
  for(int i = 0; i < num_triangles; i++)
  {

  }
}

void sphere_intersection(struct Vector3 direction, struct Color& color)
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

    if(delta >= 0)
    {
      color.r = 0;
      color.g = 0;
      color.b = 0;
    }

    float t0; 
    float t1;
  }
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
