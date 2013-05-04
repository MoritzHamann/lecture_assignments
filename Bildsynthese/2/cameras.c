/*
* @author Lena Gieseke, Alexandros Panagiotidis
* Edited by Moritz Hamann, Boitumelo Ruf
*/



#version 330
#define NUM_SPHERES (5)

layout(location = 0) out vec4 frag_color;
uniform vec2 resolution;
uniform float dof_aperture_number;



// data structures
struct Plane
{
	vec3 p0;
	vec3 r1;
	vec3 r2;
	vec4 color;
};

struct Sphere
{
	vec3 center;
	float radius;
	vec4 color;
};

struct Ray
{
	vec3 origin;
	vec3 direction;
};

// global variables
Plane g_plane;
Sphere g_spheres[NUM_SPHERES];


vec3 g_eye = vec3(0.0, 0.0, 1.0);
vec3 g_forward;
vec3 g_right;
vec3 g_up;
vec3 g_lookat;


// build world functions
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

Ray makeRay(vec3 origin, vec3 direction)
{
	Ray ray;
	
	ray.origin = origin;
	ray.direction = normalize(direction);
	
	return ray;
}

// init world functions
void initWorld()
{
	g_plane = makePlane(vec3(0.0, -2.0, -10.0)
						, vec3(1.0, 0.0, 0.0)
						, vec3(0.0, 0.0, 1.0)
						, vec4(1.0, 1.0, 0.0, 1.0)
						);
	g_spheres[0] = makeSphere(vec3(-0.5, 0.0, -1.0)
								, 0.2
								, vec4(1.0, 0.0, 0.0, 1.0)
								);

	g_spheres[1] = makeSphere(vec3(0.0, 0.0, -2.0)
								, 0.5
								, vec4(0.0, 1.0, 0.0, 1.0)
								);

	g_spheres[2] = makeSphere(vec3(0.5, 0.0, -1.5)
								, 0.35
								, vec4(0.0, 0.0, 1.0, 1.0)
								);

	g_spheres[3] = makeSphere(vec3(-1.0, g_plane.p0.y, -2.0)
								, 0.3
								, vec4(1.0, 1.0, 0.0, 1.0)
								);

	g_spheres[4] = makeSphere(vec3(1.0, 0.0, -1.0)
								, 0.3
								, vec4(1.0, 0.0, 1.0, 1.0)
								);
}

void initCamera()
{


	// Default
	g_forward = vec3(0.0, 0.0, -1.0);
	g_right = vec3(1.0, 0.0, 0.0) * float(resolution.x) / float(resolution.y);
	g_up = vec3(0.0, 1.0, 0.0);
	
	
	
	g_lookat = vec3(0.0, 0.0, 0.0);
	float fov = 45.0 / 180.0 * 3.1415962;
	float ratio = float(resolution.x) / float(resolution.y);
	vec3 forward_up_plane = vec3(0.0, 1.0, 0.0);

 	// Aufgabe 2.1 (2. Teilaufgabe): Berechnen Sie g_forward, g_right und g_up aus g_lookat, fov und ratio
	
	float tan_fov = tan(fov);
    
	g_forward = normalize(g_lookat - g_eye);
	g_right = normalize(cross(g_forward,forward_up_plane))*tan_fov*ratio;
	g_up = normalize(cross(g_forward, g_right))*tan_fov;
	
	// Ende Aufgabe 2.1 (2. Teilaufgabe)

}


// intersection functions

bool intersectPlane(Ray ray, Plane plane, out float distance)
{
	vec3 normal = normalize(cross(plane.r1, plane.r2));

	float z = dot(ray.origin - plane.p0, normal);
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

// ray tracing functions
Ray getRayOrtho(float x, float y)
{
	return makeRay(g_eye + (((2.0 * x / resolution.x) - 1.0) * g_right) + 
					(((2.0 * y / resolution.y) - 1.0) * g_up)
					, g_forward
					);
}


Ray getRayPerspective(float x, float y)
{
	// Aufgabe 2.1 (1. Teilaufgabe): Perspektivische Strahlerzeugung. Rufen Sie getRayPerspective in getColor auf anstatt getRayOrtho.
	
	return makeRay(g_eye
					, (g_forward + ( 2.0*((x+0.5)/resolution.x) - 1.0) *
					g_right - (2.0*((y+0.5)/resolution.y) - 1.0) * g_up )
					);
					
	// Ende Aufgabe 2.1 (1. Teilaufgabe)
}

vec4 radiance(Ray ray, vec3 eye)
{
  vec4 color = vec4(0.0, 0.0, 0.0, 1.0);

	float min_distance = 1000000000.0;

	float plane_distance = 0.0;
	
	if (intersectPlane(ray, g_plane, plane_distance))
	{
		min_distance = plane_distance;
		vec3 p = ray.origin + min_distance * ray.direction;

		float d1 = dot(g_plane.p0, g_plane.r1);
		float d2 = dot(g_plane.p0, g_plane.r2);

		float u = dot(p, g_plane.r1) - d1;
		float v = dot(p, g_plane.r2) - d2;

		float ui = mod(mod(u, 1.0) + 1.0, 1.0);
		float vi = mod(mod(v, 1.0) + 1.0, 1.0);

		if (ui > 0.5 ^^ vi > 0.5)
		{
			color = vec4(0.8, 0.8, 0.8, 1.0);
		}
		else
		{
			color = vec4(0.0, 0.0, 0.0, 1.0);
		}
	}
	else
	{
		color = vec4(0.0, 0.0, 0.2, 1.0);
	}

	for (int i = 0; i < NUM_SPHERES; i++)
	{
		float t1, t2;
		
		if (intersectSphere(ray, g_spheres[i], t1, t2))
		{
			float distance = min(t1, t2);
			
			if (distance > 0.0 && distance < min_distance)
			{
				min_distance = distance;
				
				vec3 p = ray.origin + min_distance * ray.direction;
		
				vec3 n = normalize(p - g_spheres[i].center);
				vec3 l = normalize(vec3(1.0, 5.0, 0.0) - p);
				vec3 v = normalize(eye - p);
				vec3 h = normalize(v + l);

				float angle = dot(n, l);
				float spec = pow(max(0.0, dot(h, n)), 50.0);

				// Ambient
				color = vec4(0.01, 0.01, 0.01, 1.0);

				// Spekular
				color += vec4(1.0) * spec;

				// Diffus
				color += g_spheres[i].color * clamp(angle, 0.0, 1.0);
			}
		}
	}

	color.a = 1.0;
  return color;
}


vec4 getColor(float x, float y)
{
	//Ray ray = getRayOrtho(x, y);	
	
	// Aufgabe 2.1
	Ray ray = getRayPerspective(x, y);
	
	vec4 color = radiance(ray, g_eye);
	
	return color;
}

// Rufen Sie diese Funktion auf anstatt getColor fuer Aufgabe 2.2
vec4 getColor_AA(float x, float y)
{
	// Aufgabe 2.2: Supersampling mit Anti-Aliasing
	
    /* Gauß Kern der Größe 7 */
    // const float stepsize = (1.0/3.0);
    // float gaussKernel[7] = float[]( 1.0,  6.0,  15.0,  20.0,  15.0,  16.0,  1.0);
	// float gaussFactor = 4096.0;
	
    /* Gauß Kern der Größe 9 */
    // const float stepsize = (1.0/4.0);
    // float gaussKernel[9] = float[]( 1.0,  8.0,  28.0,  56.0,  70.0,  56.0,  28.0, 8.0, 1.0);
	// float gaussFactor = 65535.0;

    /* Gauß Kern der Größe 11 */
    const float stepsize = (1.0/5.0);
    float gaussKernel[11] = float[]( 1.0,  10.0,  45.0,  120.0,  210.0,  252.0,  210.0, 120.0, 45.0, 10, 1.0);
	float gaussFactor = 1048576.0;
    
	vec4 color = vec4(0.0, 0.0, 0.0, 1.0);
	
	int m = 0;
	for( float i = -1; i <= 1; i = i+stepsize )
	{
		int n = 0;
		for( float j = -3.0*stepsize; j <= 3.0*stepsize; j = j+stepsize )
		{
			Ray ray = getRayPerspective((x+i), (y+j));
		
			color += (radiance(ray, g_eye) * 
				(gaussKernel[m]*gaussKernel[n])/gaussFactor);
			
			n++;
		}
		m++;
	}

	// Ende Aufgabe 2.2
	return color;
}


// Rufen Sie diese Funktion auf anstatt getColor fuer Aufgabe 2.3
vec4 getColor_DoF(float x, float y)
{
  // Aufgabe 2.3: Depth of Field
  vec4 color = vec4(0.0, 0.0, 0.0, 1.0);
  
  const float pi = 3.14159265359;
  const float focal_length = 0.6;
  const float f = 1.0;
  const float r_step = 0.01;
  const float omega_step = 2*pi/8;

  
  // Aperture opening and distance of E2
  float D = f/dof_aperture_number;
  float E2_dist = 1/(1/focal_length - 1/f);
  
  // use maximum for D (for performance)
  if (D > 3) {
    D = 3;
  }
  
  // ray through middle of lens
  Ray base_ray = getRayPerspective(x,y);
  vec3 intersection = vec3(0.0, 0.0, 0.0);
  
  // make second focusation plane of lens 
  Plane E2 = makePlane(g_eye+E2_dist*normalize(g_forward), normalize(g_up), normalize(g_right), vec4(0.0, 0.0, 0.0, 0.0));
  
  // intesection of base_ray and E2
  float distance = 0.0;
  
  int i = 1;
  if(intersectPlane(base_ray, E2, distance)){
        
    intersection = base_ray.origin + base_ray.direction * distance;
    color = radiance(base_ray, g_eye);
    

    for (float omega = 0; omega < 2*pi; omega = omega + omega_step){
      for (float r = r_step; r < D/2; r = r + r_step){
        float x_offset = r * sin(omega);
        float y_offset = r * cos(omega);
        vec3 origin = g_eye + x_offset * normalize(g_right) + y_offset * normalize(g_up);
        Ray temp = makeRay(origin, intersection - origin);
        color = color + radiance(temp, origin);
        i = i+1;
      }
    }
  }
    // Ende Aufgabe 2.3
  color = color/i;
  color.a = 1;
  return color;
}






vec4 gamma(vec4 color)
{
        // Gamma-Korrektur
        clamp(color.rgb, 0.0, 1.0);

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




void main(void)
{

	initWorld();
	initCamera();
	
	//Aufgabe 2.1
	//vec4 color = getColor(gl_FragCoord.x + 0.5, gl_FragCoord.y + 0.5);
	//Aufgabe 2.2
	//vec4 color = getColor_AA(gl_FragCoord.x + 0.5, gl_FragCoord.y + 0.5);
	//Aufgabe 2.3
	vec4 color = getColor_DoF(gl_FragCoord.x + 0.5, gl_FragCoord.y + 0.5);
	
	
	color = gamma(color);
	frag_color = color;
}
