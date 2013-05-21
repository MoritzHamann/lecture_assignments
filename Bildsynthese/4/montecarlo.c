/*
 * @author Sebastian Koch, Lena Gieseke, Alexandros Panagiotidis
 */

#version 330
#extension GL_EXT_gpu_shader4 : enable

// precision highp float;

layout(location = 0) out vec4 frag_color;

uniform vec2 resolution;
uniform sampler2D envMapTexture;
uniform vec2 lookAtAngles;
uniform int noiseSeed = 13;

uniform int sphere_switch;
uniform int BRDF_number;

uniform int samples;
uniform float shininess;
uniform vec3 color_specular;
uniform vec3 color_spheres;

// Eye of the Beholder
vec3 g_eye = vec3(0.0, 0.0, 0.0);
vec3 g_forward;
vec3 g_right;
vec3 g_up;
vec3 g_lookat;

/////////////////////////////// Random Number Generator
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
unsigned int TausStep(inout unsigned int z, 
                      int S1,
                      int S2,
                      int S3,
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
            *float( // Periods
                    TausStep(z1, 13, 19, 12, 4294967294U) ^  // p1=2^31-1
                    TausStep(z2, 2, 25, 4, 4294967288U) ^    // p2=2^30-1
                    TausStep(z3, 3, 11, 17, 4294967280U) ^   // p3=2^28-1
                    LCGStep(z4, 1664525UL, 1013904223U)        // p4=2^32
                    );
}

/////////////////////////////// End - Random Number Generator 

#define NUM_SPHERES (10)

#define M_PI ( 2.0 * asin(1.0))

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


Plane g_plane;
Sphere g_spheres[NUM_SPHERES];

Ray makeRay(vec3 origin, vec3 direction)
{
    Ray ray;

    ray.origin = origin;
    ray.direction = normalize(direction);

    return ray;
}


Plane makePlane(vec3 p0, vec3 r1, vec3 r2, vec4 color)
{
    Plane g_plane;

    g_plane.p0 = p0;
    g_plane.r1 = normalize(r1);
    g_plane.r2 = normalize(r2);
    g_plane.color = color;

    return g_plane;
}


Sphere makeSphere(in vec3 center, in float radius, in vec4 color)
{
    Sphere sphere;

    sphere.center = center;
    sphere.radius = radius;
    sphere.color = color;

    return sphere;
}


bool intersectPlane(Ray ray, Plane g_plane, out float distance)
{
    vec3 normal = normalize(cross(g_plane.r1, g_plane.r2));

    float z = dot(ray.origin - g_plane.p0, normal);
    float n = dot(ray.direction, normal);

    distance = 0.0;

    if (abs(n) > 1.0e-6)
    {
        distance = -(z / n);
        return distance >= 0.0;
    }

    return false;
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

    if (det < 0.0)
    {
        return false;
    }

    det = sqrt(det);

    float q;

    if (b < 0.0)
    {
        q = -0.5 * (b - det);
    }
    else
    {
        q = -0.5 * (b + det);
    }

    t1 = q / a;
    t2 = c / q;

    return true;
}

vec3 sphericalToCartesian(vec3 spherical) 
{

    float phi = spherical.x;
    float theta = spherical.y;
    float r = spherical.z;

    vec3 cart;

    cart.x = sin(theta) * cos(phi);
    cart.y = cos(theta);
    cart.z = sin(theta) * sin(phi);

    return cart;
}


vec3 cartesianToSpherical(vec3 cart) 
{

    float phi = atan(cart.x, cart.z);
    float theta = acos(cart.y);
    float r = length(cart);

    return vec3(phi, theta, r);
}

vec3 sampleEnvMap(vec3 dir)
{
    vec3 color = vec3( 0.0);

    vec3 angles;
    vec2 coords;

    angles = cartesianToSpherical( dir);

    color = texture2D( envMapTexture,
                       vec2( 1.0 - angles.x/( 2.0 * M_PI), angles.y / M_PI)).xyz;

    return color;
}

vec3 compProd(vec3 a, vec3 b) 
{
    return vec3(a.r * b.r, a.g * b.g, a.b * b.b);
}


/////////////////////////////// BRDFs

// n = normal; v = ray.direction; l = light.direction
vec3 phongBRDF( vec3 n,
                vec3 v,
                vec3 l,
                vec3 diffuseCol,
                vec3 specularCol,
                float exponent)
{

    vec3 r = reflect(-v, n);
	
    // diffuse term
    vec3 diff = diffuseCol / M_PI;

    // specular term
    float specFac = 10.0;
    specFac *= pow(max(0.0, dot(r, l)), exponent);

    return diff + specularCol * specFac;
}

vec3 blinnPhongBRDF( vec3 n,
                     vec3 v,
                     vec3 l,
                     vec3 diffuseCol,
                     vec3 specularCol,
                     float exponent)
{

    vec3 h = normalize(l + v);

    // diffuse term
    vec3 diff = diffuseCol / M_PI;

    // specular term
    float specFac = (exponent + 2.0) / ( 2.0 * M_PI);
    specFac *= pow( dot(h, n), exponent);

    return diff + specularCol * specFac;
}


vec3 lambertBRDF(vec3 diffuse){
    return diffuse/M_PI;
}




/////////////////////////////// Uniform Sampling

vec3 getUniformRandomSampleDirectionUpper(in vec3 N){

    // TODO 4.1.1
	// Uniform Sampling of the hemisphere
	
    vec3 sampleDirection = vec3(0.0);

    float z = HybridTaus();	
    float phi = (2*M_PI)*z;
    
    sampleDirection.x = sqrt(1-pow(z,2))*sin(phi);
    sampleDirection.y = sqrt(1-pow(z,2))*cos(phi);
    sampleDirection.z = z;
    
    if(dot(sampleDirection, N) < 0)
    {
        sampleDirection = -sampleDirection;
    }
	
    return sampleDirection;
	// END 4.1.1
}



vec3 uniformSampleLambert(in vec3 N, in vec3 diffuse)
{
    // TODO 4.1.2 : 
	// Uniform Sampling, Lambert BRDF
    vec3 color = color_spheres;
    vec3 sampleDirection = vec3(0.0);
    
    for(int i = 0; i < samples; i++)
    {
        sampleDirection = getUniformRandomSampleDirectionUpper(N);
        color += (sampleEnvMap(sampleDirection) * lambertBRDF(diffuse));
    }
	color /= float(samples);
    
	return color;
	// END 4.1.2
}

vec3 uniformSamplePhong(in vec3 V, in vec3 N, in vec3 diffuse,
                        in vec3 specular, in float shiny)
{

    // TODO 4.1.3
    ////Uniform Sampling, Phong BRDF
    vec3 color = color_spheres;
    vec3 sampleDirection = vec3(0.0);
    
    for(int i = 0; i < samples; i++)
    {
        sampleDirection = getUniformRandomSampleDirectionUpper(N);
        color += (sampleEnvMap(sampleDirection) * phongBRDF(N, V, sampleDirection, diffuse, specular, shiny));
    }
	color /= float(samples);
    
	return color;
	// END 4.1.3
}



/////////////////////////////// IMPORTANCE SAMPLING

vec3 orthogonalVector(vec3 n) {

    if((abs(n.y) >= 0.9*abs(n.x)) && (abs(n.z) >= 0.9*abs(n.x)))
        return vec3(0.0, -n.z, n.y);
    else
        if((abs(n.x) >= 0.9*abs(n.y)) && (abs(n.z) >= 0.9*abs(n.y)))
            return vec3(-n.z, 0.0, n.x);
        else
            return vec3(-n.y, n.x, 0.0);
}


void localFrame(in vec3 n, out vec3 a, out vec3 b) {

    a = normalize(orthogonalVector(n));
    b = normalize(cross(a, n));
}


mat3 localTransformation(vec3 n) {

    vec3 a, b;
    localFrame(n, a, b);

    mat3 A = mat3(
                b.x, b.y, b.z, // first column (not row!)
                a.x, a.y, a.z, // second column
                n.x, n.y, n.z  // third column
                );

    return A;
}

vec3 cosineHemisphereSample3(float u1, float u2, float e) 
{

    // TODO 4.2.1
	// Vector constuction based on the cosinus distribution
	vec3 dir = vec3(0.0);
    
    float exp = 2/(e+1);

    dir.x = cos(2*M_PI*u1)*sqrt(1-pow(u2,(2/(e+1))));
    dir.y = sin(2*M_PI*u1)*sqrt(1-pow(u2,(2/(e+1))));
	dir.z = pow(u2,(1/(e+1)));
    
	return dir;
	// END 4.2.1
}
    
vec3 cosineHemisphereSample1(float e) 
{
    vec3 dir = vec3(0.0);

    float u1 = HybridTaus();
    float u2 = HybridTaus();
    
    dir = cosineHemisphereSample3(u1, u2, e);
    
    return dir;
}

float cosineHemispherePDF(float cos_theta, float e) 
{

    // TODO 4.2.2
    // Pdf calculation for cosinus distribution
	float pdf = 0.0;
	
    pdf = (e+1)/(2*M_PI) * pow(cos_theta,e);
	
	return pdf;
	// END 4.2.2
}


vec3 importanceSampleLambert(in vec3 N, in vec3 diffuse)
{

    // TODO 4.2.3
    // Hemisphere with cos^e importance sampling, Lambert BRDF
	vec3 color = color_spheres;
    vec3 sampleDirection = vec3(0.0);
    vec3 transSampleDirection = vec3(0.0);
    float exp = 1.0;
    float cos_theta = 0.0;
    
    mat3 transMatrix = localTransformation(N);
	
    for(int i = 0; i < samples; i++)
    {
        sampleDirection = cosineHemisphereSample1(exp);
        transSampleDirection = transMatrix * sampleDirection;
        cos_theta = dot(N,sampleDirection)/(length(N)*length(sampleDirection));
        
        color += (sampleEnvMap(sampleDirection) * lambertBRDF(diffuse))/cosineHemispherePDF(cos_theta, exp);
        // vec3 debug = vec3(cosineHemispherePDF(cosTheta, exp));
        // color = debug;
        // color /= cosineHemispherePDF(cos_theta, exp);
	}
    
	return color;
	// END 4.2.3
}




vec3 importanceSamplePhong(in vec3 V, in vec3 N, in vec3 diffuse,
                           in vec3 specular, in float shiny)
{

    // ***************************BONUS*******************************
    // DO YOU DARE IT??

    // TODO Aufgabe 4.2.4
    // Hemisphere with cos^e importance sampling, Phong BRDF
	
	vec3 color = color_spheres;

    // ...
	
	return color;
}


vec4 radiance(Ray ray, vec3 eye)
{


    vec4 color = vec4(0.0, 0.0, 0.0, 0.0);

    float minDistance = 1000000000.0;

    float planeDistance = 0.0;

    int hit_sphere = 0;
    // ray-sphere-intersection
    for (int i = 0; i < NUM_SPHERES; i++){
        float t1, t2;

        if (intersectSphere(ray, g_spheres[i], t1, t2)){

            float distance = min(t1, t2);

            if (distance > 0.0 && distance < minDistance){
                hit_sphere = 1;
                minDistance = distance;

                vec3 p = ray.origin + minDistance * ray.direction;

                vec3 n = normalize(p - g_spheres[i].center);

                vec3 l = normalize(vec3(1.0, 5.0, 0.0) - p);
                vec3 v = normalize(eye - p);


                for(int j = 0; j < samples; ++j){


                    switch(BRDF_number){
                    case(0):
                        color.rgb += uniformSampleLambert(n, color_spheres);
                        break;
                    case(1):
                        color.rgb += importanceSampleLambert(n, color_spheres);
                        break;
                    case(2):
                        color.rgb += uniformSamplePhong(v, n,
                                                        color_spheres,
                                                        color_specular, shininess);
                        break;
                    case(3):
                        color.rgb += importanceSamplePhong(v, n,
                                                           color_spheres,
                                                           color_specular, shininess);
                        break;
                    default:
                        color += g_spheres[i].color;
                        //color += vec4(1.0,0.0,0.0,1.0);
                    }
                }

                color /= float( samples);

            }
        }
    }

    if(hit_sphere == 0)
	{
        color = vec4( sampleEnvMap( ray.direction), 1.0);
    }

    color.a = 1.0;
    return color;
}


vec4 gamma(vec4 color)
{
    // Gamma-Correction

    clamp(color.rgb, 0.0, 1.0);

    const float magic = 0.0031308;
    const float factor = 12.92;
    const float a = 0.055;
    const float aPlus1 = a + 1.0;
    const float inv24 = 1.0 / 2.4;

    color.r = (color.r <= magic) ?
                (color.r * factor) : (aPlus1 * pow(color.r, inv24) - a);
    color.g = (color.g <= magic) ?
                (color.g * factor) : (aPlus1 * pow(color.g, inv24) - a);
    color.b = (color.b <= magic) ?
                (color.b * factor) : (aPlus1 * pow(color.b, inv24) - a);
    color.a = 1.0;

    return color;
}


void initWorld()
{
    g_plane = makePlane(vec3(0.0, -2.0, -10.0)
                        , vec3(1.0, 0.0, 0.0)
                        , vec3(0.0, 0.0, 1.0)
                        , vec4(1.0, 1.0, 0.0, 1.0)
                        );

    vec4 spcol[NUM_SPHERES];
    spcol[0] = vec4(1.0, 1.0, 1.0, 1.0);
    spcol[1] = vec4(1.0, 1.0, 0.0, 1.0);
    spcol[2] = vec4(1.0, 0.0, 1.0, 1.0);
    spcol[3] = vec4(0.0, 0.0, 1.0, 1.0);
    spcol[4] = vec4(0.0, 1.0, 0.5, 1.0);
    spcol[5] = vec4(0.5, 1.0, 0.3, 1.0);
    spcol[6] = vec4(1.0, 0.5, 0.2, 1.0);
    spcol[7] = vec4(0.3, 1.0, 0.4, 1.0);
    spcol[8] = vec4(0.4, 0.2, 0.8, 1.0);
    spcol[9] = vec4(0.0, 1.0, 0.0, 1.0);


    g_spheres[0] = makeSphere(vec3(0.0, 0.0, -1.0) * 2.0
                              , 0.5
                              , spcol[0]
            );

    g_spheres[1] = makeSphere(vec3(-1.0, 0.0, 0.0) * 2.0
                              , 0.5
                              , spcol[1]
            );

    g_spheres[2] = makeSphere(vec3(1.0, 0.0, 0.0) * 2.0
                              , 0.5
                              , spcol[2]
            );

    g_spheres[3] = makeSphere(vec3(0.0, 1.0, 0.0) * 2.0
                              , 0.5
                              , spcol[3]
            );

    g_spheres[4] = makeSphere(vec3(0.0, -1.0, 0.0) * 2.0
                              , 0.5
                              , spcol[4]
            );

    g_spheres[5] = makeSphere(vec3(0.0, 0.0, 1.0) * 2.0
                              , 0.5
                              , spcol[5]
            );

    g_spheres[6] = makeSphere(normalize(vec3(1.0, 0.5, -1.0)) * 2.0
                              , 0.5
                              , spcol[6]
            );

    g_spheres[7] = makeSphere(normalize(vec3(1.0, 1.0, 1.0)) * 2.0
                              , 0.5
                              , spcol[7]
            );

    g_spheres[8] = makeSphere(normalize(vec3(-1.0, -0.5, 1.0)) * 2.0
                              , 0.5
                              , spcol[8]
            );

    g_spheres[9] = makeSphere(normalize(vec3(-1.0, -1.0, -1.0)) * 2.0
                              , 0.5
                              , spcol[9]
            );

    //helpful for debugging: all spheres' colors to white
    //for (int i = 0; i < NUM_SPHERES; i++){
    	//g_spheres[i].color = vec4(1.0);
    //}


}


void initCamera()
{
    vec3 vert = vec3(0.0, 1.0, 0.0);

    g_lookat = sphericalToCartesian( vec3(lookAtAngles, 1.0) );

    float FoV = 45.0 / 180.0 * 3.1415962;
    float ratio = float(resolution.x) / float(resolution.y);

    // Default
    g_forward = vec3(0.0, 0.0, -1.0);
    g_right = vec3(1.0, 0.0, 0.0) * float(resolution.x) / float(resolution.y);
    g_up = vec3(0.0, 1.0, 0.0);


    float tanFoV = 2.0 * tan(FoV * 0.5);

    g_forward = normalize( g_lookat);
    g_right = normalize(cross(g_forward, vert)) * tanFoV * ratio;
    g_up = normalize(cross(g_right, g_forward)) * tanFoV;


}


Ray getRayOrtho(float x, float y)
{
    return makeRay(g_eye + (((2.0 * x / resolution.x) - 1.0) * g_right)
                   + (((2.0 * y / resolution.y) - 1.0) * g_up)
                   , g_forward
                   );
}


Ray getRayPerspective(float x, float y)
{

    return makeRay(g_eye
                   , g_forward + (((2.0 * x / resolution.x) - 1.0) * g_right)
                   + (((2.0 * y / resolution.y) - 1.0) * g_up)
                   );
}


vec4 getColor(float x, float y)
{
    Ray ray = getRayPerspective(x, y);
    vec4 color = radiance(ray, g_eye);

    return color;
}

void main(void)
{
    initWorld();
    initCamera();
    initRNG();

    vec4 color = getColor(gl_FragCoord.x, gl_FragCoord.y);
    color = gamma(color);

    frag_color = color;
}
