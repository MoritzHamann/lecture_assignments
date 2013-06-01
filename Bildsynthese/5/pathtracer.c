/*
 * @author Sebastian Koch, Lena Gieseke, Alexandros Panagiotidis
 */

#version 330
#extension GL_EXT_gpu_shader4 : enable

layout(location = 0) out vec4 frag_color;

uniform vec2 resolution;
uniform sampler2D envMapTexture;
uniform vec2 lookAtAngles;

uniform int num_rays;
uniform int num_bounces;
uniform int path_tracing;



// Eye of the Beholder
vec3 I = vec3(0.0, 0.0, 0.0);
vec3 forward;
vec3 right;
vec3 up;
vec3 lookAt;

float distanceToCamera = 0.0;

const int numSamples = 64;

#define NUM_SPHERES 2
#define NUM_TRIANGLES 8

#define M_PI 3.1415926535897932384626433832795


struct Ray
{
  vec3 origin;
  vec3 direction;
};


struct Sphere
{
  vec3 center;
  float radius;
  vec4 color;
};


struct Plane
{
  vec3 p0;
  vec3 r1;
  vec3 r2;
  vec4 color;
};


struct Triangle
{
  vec3 v0;
  vec3 v1;
  vec3 v2;
  vec3 normal;
  vec4 color;
};


Sphere spheres[NUM_SPHERES];
Triangle triangles[NUM_TRIANGLES];


Ray makeRay(vec3 origin, vec3 direction)
{
  Ray ray;

  ray.origin = origin;
  ray.direction = normalize(direction);

  return ray;
}


Plane makePlane(vec3 p0, vec3 r1, vec3 r2, vec4 color)
{
  Plane plane;

  plane.p0 = p0;
  plane.r1 = normalize(r1);
  plane.r2 = normalize(r2);
  plane.color = color;

  return plane;
}


Sphere makeSphere(in vec3 center, in float radius, in vec4 color)
{
  Sphere sphere;

  sphere.center = center;
  sphere.radius = radius;
  sphere.color = color;

  return sphere;
}


Sphere translateSphere(Sphere sphere, vec3 t)
{
  vec3 center = sphere.center + t;

  return makeSphere( center, sphere.radius, sphere.color);
}


Triangle makeTriangle(in vec3 v0, in vec3 v1, in vec3 v2, in vec4 color)
{
  Triangle triangle;

  triangle.v0 = v0;
  triangle.v1 = v1;
  triangle.v2 = v2;
  triangle.normal = normalize(cross(v0 - v2, v1 - v2));
  triangle.color = color;

  return triangle;
}


Triangle scaleTriangle(Triangle triangle, float s)
{
  vec3 normal = triangle.normal;

  vec3 v0 = triangle.v0 * s;
  vec3 v1 = triangle.v1 * s;
  vec3 v2 = triangle.v2 * s;

  Triangle tria = makeTriangle(v0, v1, v2, triangle.color);
  tria.normal = triangle.normal;
  return tria;
}


Triangle translateTriangle(Triangle triangle, vec3 t)
{
  vec3 normal = triangle.normal;

  vec3 v0 = triangle.v0 + t;
  vec3 v1 = triangle.v1 + t;
  vec3 v2 = triangle.v2 + t;

  Triangle tria = makeTriangle(v0, v1, v2, triangle.color);
  tria.normal = triangle.normal;
  return tria;
}


bool intersectPlane(Ray ray, Plane plane, out float distance)
{
  vec3 normal = normalize(cross(plane.r1, plane.r2));

  float z = dot(ray.origin - plane.p0, normal);
  float n = dot(ray.direction, normal);

  distance = 0.0;

  if (abs(n) > 1.0e-6){
      distance = -(z / n);
      return distance >= 0.0;
  }

  return false;
}

//see: http://www.blackpawn.com/texts/pointinpoly/default.html
bool pointOnSameSide(vec3 p1, vec3 p2, vec3 a, vec3 b)
{
    vec3 cp1 = cross((b-a), (p1-a));
    vec3 cp2 = cross((b-a), (p2-a));
    
    if(dot(cp1, cp2)>= 0)
    {
        return true;
    }
    else
    {
        return false;
    }
}

//see: http://www.blackpawn.com/texts/pointinpoly/default.html
bool pointInTriangle(vec3 p, vec3 a, vec3 b, vec3 c)
{
    if(pointOnSameSide(p,a,b,c) && pointOnSameSide(p,b,a,c) && pointOnSameSide(p,c,a,b))
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool intersectTriangle(Ray ray, Triangle triangle, out float distance)
{
  bool retVal = false;
  float tmpDistance = 0.0;
  distance = 0.0;
  
  // erstellen der entsprechenden Ebene zum Dreieck
  Plane plane = makePlane(triangle.v2, (triangle.v0-triangle.v2), (triangle.v1-triangle.v2), triangle.color);
  
  // prüfen, ob Strahl die Ebene schneidet.
  if(intersectPlane(ray,plane,tmpDistance))
  {
    // wird die Ebene vom strahl geschnitten prüfe, ob der schnittpunkt innerhalb des Dreiecks liegt.
    vec3 intersectionPoint = ray.origin + (tmpDistance*ray.direction);
    
    // Um zu überprüfen ob ein Punkt in einem Dreieck liegt verfolgen wir den Ansatz auf: http://www.blackpawn.com/texts/pointinpoly/default.html
    if(pointInTriangle(intersectionPoint, triangle.v0, triangle.v1, triangle.v2))
    {
        retVal = true;
        distance = tmpDistance;
    }    
  }
  
  // wird die Ebene vom Strahl nicht geschnitten, wird das Dreieck auch nicht geschnitten.
  
  return retVal;

}


bool intersectSphere(in Ray ray, in Sphere sphere, out float t1, out float t2)
{
  vec3 p = ray.origin - sphere.center;

  float a = dot(ray.direction, ray.direction);
  float b = 2.0 * dot(p, ray.direction);
  float c = dot(p, p) - (sphere.radius * sphere.radius);

  float det = b * b - 4.0 * a * c;

  t1 = 0.0;
  t2 = 0.0;

  if (det < 0.0){
    return false;
  }

  det = sqrt(det);

  float q;

  if (b < 0.0){
      q = -0.5 * (b - det);
  }
  else{
      q = -0.5 * (b + det);
    }

  t1 = q / a;
  t2 = c / q;

  return true;
}


vec3 sphericalToCartesian(vec3 spherical) {

  float phi = spherical.x;
  float theta = spherical.y;
  float r = spherical.z;

  vec3 cart;

  cart.x = sin(theta) * sin(phi);
  cart.y = cos(theta);
  cart.z = -sin(theta) * cos(phi);

  return r * cart;
}


vec3 cartesianToSpherical(vec3 cart) {

  float phi = atan(cart.x, -cart.z);
  float theta = acos(cart.y);
  float r = length(cart);

  return vec3(phi, theta, r);
}


vec3 sampleEnvMap(vec3 dir)
{
  vec3 color = vec3( 0.5);
  vec3 angles;

  angles = cartesianToSpherical( dir);
  angles.x += lookAtAngles.x;
  color = texture2D( envMapTexture,
                     vec2( angles.x/( 2.0 * M_PI), angles.y / M_PI)).rgb;

  return color;
}


vec3 compProd(vec3 a, vec3 b) {
  return vec3(a.r * b.r, a.g * b.g, a.b * b.b);
}

/////////////////////// Random Number Generator ///////////////////////
unsigned int z1, z2, z3, z4;

void initRNG()
{
        z1 = 128U + uint(gl_FragCoord.x);
        z2 = 128U + uint(gl_FragCoord.y);
        z3 = 128U + uint(gl_FragCoord.x + gl_FragCoord.y);
        z4 = 128U + uint(gl_FragCoord.x * gl_FragCoord.y);
}

// S1, S2, S3, and M are all constants, and z is part of the
// private per-thread generator state.
unsigned int TausStep(inout unsigned int z, int S1, int S2, int S3,
                      unsigned int M) {
  unsigned int b=(((z << uint(S1)) ^ z) >> uint(S2));
  return z = (((z & M) << uint(S3)) ^ b);
}

// A and C are constants
unsigned int LCGStep(inout unsigned int z, unsigned int A, unsigned int C) {
  return z=(A*z+C);
}

float HybridTaus() {
  // Combined period is lcm(p1,p2,p3,p4)~ 2^121
  return 2.3283064365387e-10
    * float(         // Periods
            TausStep(z1, 13, 19, 12, 4294967294) ^  // p1=2^31-1
            TausStep(z2, 2, 25, 4, 4294967288) ^    // p2=2^30-1
            TausStep(z3, 3, 11, 17, 4294967280) ^   // p3=2^28-1
            LCGStep(z4, 1664525U, 1013904223U)        // p4=2^32
                     );
}



bool intersectScene(Ray ray, out vec4 color, out vec3 p, out vec3 n) {
  float minDistance = 1000.0;
  color = vec4(0.0);
  p = vec3(0.0);
  n = vec3(0.0);
  // ray-triangle-intersection
  for (int i = 0; i < NUM_TRIANGLES; i++) {
    float d;

    if (intersectTriangle(ray, triangles[i], d) &&
            d > 0.0 && d < minDistance) {
      color = triangles[i].color;
      minDistance = d;
      p = ray.origin + minDistance * ray.direction;
      n = triangles[i].normal;
    }
  }

  // ray-sphere-intersection
  for (int i = 0; i < NUM_SPHERES; i++) {
    float t1, t2;
    if (intersectSphere(ray, spheres[i], t1, t2)) {
      float distance = min(t1, t2);
      if (distance > 0.0 && distance < minDistance) {
        color = spheres[i].color;
        minDistance = distance;
        p = ray.origin + minDistance * ray.direction;
        n = normalize(p - spheres[i].center);
      }
    }
  }

  // we hit something!
  if( minDistance < 999.9 )  {
    return true;
  }
  else
    return false;
}


vec3 uniformSampleDirection(void) {

  float z = (HybridTaus() * 2.0) - 1.0;
  float phi = HybridTaus() * 2.0 * M_PI;
  vec3 sampleDirection = vec3(sqrt( 1.0 - z*z) * sin( phi)
                              , sqrt( 1.0 - z*z) * cos( phi)
                              , z );

  return normalize(sampleDirection);
}

vec3 lambertBRDF(vec3 diffuse){
    return diffuse/M_PI;
}

vec4 radiance(Ray ray, vec3 eye) {
  vec4 color = vec4( 0.0);
  
  // Anfang Aufgabe 5.3
  vec3 intersectionPoint = vec3(0.0);
  vec3 normal = vec3(0.0);
  vec3 newSampleDirection = vec3(0.0);
  vec4 tmpColor = vec4(1.0);
  vec3 diffuse = vec3(1.0);
  
  int i = 1;
  while(i <= num_bounces)
  {
    // Prüfe ob der Strahl die Szene Schneidet.
    if(intersectScene(ray, tmpColor, intersectionPoint, normal))
    {
        color.rgb += (tmpColor.rgb * lambertBRDF(diffuse));
            
        diffuse = tmpColor.rgb;
        
        // generiere ein Sample, das in Richtung der oberen Halbkugel zeigt.
        do
        {
            newSampleDirection = uniformSampleDirection();
        }while(dot(normal,newSampleDirection) < 0);
        
        ray.origin = intersectionPoint;
        ray.direction = newSampleDirection;
        i++;
    }
    else
    {
        if(i == 1)
            color.rgb += sampleEnvMap(ray.direction);
        else
            color.rgb += (sampleEnvMap(ray.direction) * lambertBRDF(diffuse)) ;
        break;
    }
  }
  
  if(i < num_bounces)
  {    
    color.rgb /= float(i);
  }
  else
  {
    // der Strahl hat die EnvMap nicht getroffen.
    color = vec4(0.0);
  }
  
  // Ende Aufgabe 5.3

  return color;

}


vec4 gamma(vec4 color)
{
  // Gamma-Korrektur

  const float magic = 0.0031308;
  const float factor = 12.92;
  const float a = 0.055;
  const float aPlus1 = a + 1.0;
  const float inv24 = 1.0 / 2.4;

  color.r = (color.r <= magic) ? (color.r * factor) : (aPlus1 * pow(color.r, inv24) - a);
  color.g = (color.g <= magic) ? (color.g * factor) : (aPlus1 * pow(color.g, inv24) - a);
  color.b = (color.b <= magic) ? (color.b * factor) : (aPlus1 * pow(color.b, inv24) - a);
  color.a = 1.0;

  return color;
}


void initWorld()
{
  // right wall (cyan)
  vec4 cyan = vec4(0.0, 0.8, 0.8, 1.0) / M_PI;
  triangles[0] = makeTriangle(vec3( 1.0,-1.0, -1.0),
                  vec3( 1.0, 1.0, -1.0),
                  vec3( 1.0,-1.0,  1.0),
                  cyan);

  triangles[1] = makeTriangle(vec3( 1.0, 1.0, -1.0),
                  vec3( 1.0, 1.0,  1.0) ,
                  vec3( 1.0,-1.0,  1.0),
                  cyan);

  // floor (gray)
  vec4 gray01 = vec4(0.8, 0.8, 0.8, 1.0) / M_PI;
  triangles[2] = makeTriangle(vec3(-1.0, -1.0, -1.0),
                  vec3( 1.0, -1.0, -1.0),
                  vec3(-1.0, -1.0,  1.0),
                  gray01);

  triangles[3] = makeTriangle(vec3( 1.0, -1.0, -1.0),
                  vec3( 1.0, -1.0,  1.0) ,
                  vec3(-1.0, -1.0,  1.0),
                  gray01);

  // left wall (red)
  vec4 red = vec4(0.8, 0.0, 0.0, 1.0) / M_PI;
  triangles[4] = makeTriangle(vec3( -1.0,-1.0, -1.0),
                  vec3( -1.0, 1.0, -1.0),
                  vec3( -1.0,-1.0,  1.0),
                  red);

  triangles[5] = makeTriangle(vec3( -1.0, 1.0, -1.0),
                  vec3( -1.0, 1.0,  1.0) ,
                  vec3( -1.0,-1.0,  1.0),
                  red);

  // front wall (gray)
  triangles[6] = makeTriangle(vec3(-1.0, -1.0, -1.0),
                  vec3( 1.0, -1.0, -1.0),
                  vec3(-1.0,  1.0, -1.0),
                  gray01);

  triangles[7] = makeTriangle(vec3( 1.0, -1.0, -1.0),
                  vec3( 1.0,  1.0, -1.0),
                  vec3(-1.0,  1.0, -1.0),
                  gray01);

  //// ceiling (gray)
  //triangles[8] = makeTriangle(vec3(-1.0, 1.0, -1.0),
  //                            vec3( 1.0, 1.0, -1.0),
  //							vec3(-1.0, 1.0,  1.0),
  //							gray01);

  //triangles[9] = makeTriangle(vec3( 1.0, 1.0, -1.0),
  //                            vec3( 1.0, 1.0,  1.0) ,
  //							vec3(-1.0, 1.0,  1.0),
  //							gray01);


  // spheres
  vec4 white = vec4(1.0, 1.0, 1.0, 1.0) / M_PI;
  spheres[0] = makeSphere(vec3(0.0,-0.5, 0.5),
              0.5,
              white);

  vec4 gray02 = vec4( 0.6, 0.6, 0.6, 1.0) / M_PI;
  spheres[1] = makeSphere(vec3(0.5, 0.0, 0.0),
              0.4,
              gray02);

  // -- translate all objects --
  vec3 trans = vec3(0.0, 0.0, -3.5);
  for( int i = 0; i < NUM_TRIANGLES; i++ ) {
    triangles[i] = translateTriangle( triangles[i], trans);
  }
  for( int i = 0; i < NUM_SPHERES; i++ ) {
    spheres[i] = translateSphere( spheres[i], trans);
  }

  // -- set correct normals --
  // right wall (cyan)
  triangles[0].normal = vec3(-1.0, 0.0, 0.0);
  triangles[1].normal = vec3(-1.0, 0.0, 0.0);
  // floor (gray)
  triangles[2].normal = vec3( 0.0, 1.0, 0.0);
  triangles[3].normal = vec3( 0.0, 1.0, 0.0);
  // left wall (red)
  triangles[4].normal = vec3( 1.0, 0.0, 0.0);
  triangles[5].normal = vec3( 1.0, 0.0, 0.0);
  // front wall (gray)
  triangles[6].normal = vec3( 0.0, 0.0, 1.0);
  triangles[7].normal = vec3( 0.0, 0.0, 1.0);
  // ceiling (gray)
  //triangles[8].normal = vec3( 0.0,-1.0, 0.0);
  //triangles[9].normal = vec3( 0.0,-1.0, 0.0);
}


void initCamera(float fov)
{
  const vec3 vert = vec3(0.0, 1.0, 0.0);

  lookAt = vec3( 0.0, 0.0,-1.0);
  //lookAt = sphericalToCartesian(vec3(lookAtAngles, 1.0));
  float ratio = float(resolution.x) / float(resolution.y);

  float tanFoV = tan((fov / 180.0 * 3.1415962) * 0.5);

  forward = normalize(lookAt);
  right = normalize(cross(forward, vert)) * tanFoV * ratio;
  up = normalize(cross(right, forward)) * tanFoV;
}


Ray getRayPerspective(float x, float y)
{
  return makeRay(I, forward + (((2.0 * x / resolution.x) - 1.0) * right)
                 + (((2.0 * y / resolution.y) - 1.0) * up));
}


vec4 getColor(float x, float y)
{

  vec4 color = vec4(0.0, 0.0, 0.0, 1.0);


  if(path_tracing == 0){
    Ray r = getRayPerspective( x,  y);
    vec3 p;
    vec3 n;
    vec4 c;

  
    if( intersectScene( r, c, p, n) ) {
      color = c;
    }else{
      color = vec4(sampleEnvMap( r.direction), 1.0);
    }

  }else{
    float offset_x = 0.0;
    float offset_y = 0.0;
    
    for(int i = 0; i < num_rays; i++)
    {
        offset_x = HybridTaus() - 0.5;
        offset_y = HybridTaus() - 0.5;
        
        Ray r = getRayPerspective( (x+offset_x),  (y+offset_y));
      
        color = radiance(r, r.origin);
    }
  }
  color.a = 1.0;

  return color;
}


void main(void)
{
  initWorld();
  initCamera(60.0);
  initRNG();

  vec4 color = vec4(0.0);
  color.a = 1.0;

  color = getColor(gl_FragCoord.x + 0.5, gl_FragCoord.y + 0.5);
  color = gamma(color);

  frag_color = color;
}
