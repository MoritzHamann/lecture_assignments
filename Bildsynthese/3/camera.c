/*
 * @author Sebastian Koch, Lena Gieseke, Alexandros Panagiotidis
 */

#version 330

// precision highp float;

layout(location = 0) out vec4 frag_color;

uniform vec2 resolution;
uniform sampler2D envMapTexture;
uniform vec2 lookAtAngles;
uniform float gain_number;
uniform int gain_switch;
uniform vec3 gain_values;
uniform int sphere_switch;
uniform int BRDF_number;


// Eye of the Beholder
vec3 g_eye = vec3(0.0, 0.0, 0.0);
vec3 g_forward;
vec3 g_right;
vec3 g_up;
vec3 g_lookat;



const int numLights = 64;

const vec3 lightDirections [64] = vec3[](
                                         vec3(0.65999,	-0.358729,	0.660096),
                                         vec3(0.183328,	0.869764,	-0.45815),
                                         vec3(-0.239132,	0.345878,	-0.907295),
                                         vec3(0.440225,	-0.0505455,	0.896464),
                                         vec3(-0.0180877,	-0.0822758,	0.996445),
                                         vec3(0.925856,	0.279891,	0.253873),
                                         vec3(-0.34909,	-0.103842,	-0.931318),
                                         vec3(0.235002,	-0.96023,	0.150776),
                                         vec3(-0.586574,	-0.807475,	0.0625734),
                                         vec3(-0.746192,	0.646104,	-0.160461),
                                         vec3(0.331124,	0.879095,	0.342854),
                                         vec3(-0.897824,	0.236481,	-0.371468),
                                         vec3(0.49793,	0.860148,	-0.110503),
                                         vec3(0.53136,	0.525304,	0.664615),
                                         vec3(0.0439909,	0.998015,	-0.0450684),
                                         vec3(-0.694336,	0.27298,	0.665868),
                                         vec3(-0.470362,	0.873236,	0.12735),
                                         vec3(-0.852479,	-0.464135,	0.240536),
                                         vec3(0.463292,	-0.242637,	-0.852343),
                                         vec3(-0.988212,	-0.131238,	-0.0788242),
                                         vec3(-0.744799,	0.575586,	0.337601),
                                         vec3(0.812092,	-0.182242,	-0.554341),
                                         vec3(-0.328797,	0.897185,	-0.294876),
                                         vec3(0.116558,	0.0573657,	-0.991526),
                                         vec3(0.0820923,	0.697463,	0.711903),
                                         vec3(0.756436,	-0.590919,	0.280391),
                                         vec3(0.227287,	-0.467987,	0.854007),
                                         vec3(-0.456962,	-0.0276479,	0.889056),
                                         vec3(0.660338,	-0.591154,	-0.463132),
                                         vec3(0.907803,	-0.401316,	-0.121817),
                                         vec3(0.990642,	0.0395601,	-0.130624),
                                         vec3(0.753718,	0.136227,	0.642924),
                                         vec3(-0.117144,	0.921614,	0.370007),
                                         vec3(0.829053,	0.256361,	-0.49694),
                                         vec3(0.560421,	0.629669,	-0.538001),
                                         vec3(0.834776,	0.530858,	-0.146076),
                                         vec3(-0.0250712,	-0.781426,	0.623495),
                                         vec3(-0.0161355,	-0.391658,	-0.919969),
                                         vec3(-0.823442,	-0.516925,	-0.233948),
                                         vec3(-0.694406,	-0.306917,	0.650847),
                                         vec3(-0.5037,	-0.769194,	-0.393227),
                                         vec3(-0.239166,	0.366913,	0.898986),
                                         vec3(0.410585,	-0.754856,	0.51148),
                                         vec3(-0.123284,	-0.761756,	-0.636026),
                                         vec3(0.304806,	-0.633937,	-0.710786),
                                         vec3(-0.198945,	-0.945761,	0.256821),
                                         vec3(-0.645087,	0.163386,	-0.746437),
                                         vec3(0.257978,	-0.903477,	-0.342312),
                                         vec3(-0.917414,	-0.0134328,	0.397707),
                                         vec3(0.218722,	0.323173,	0.920717),
                                         vec3(0.534413,	0.20025,	-0.821159),
                                         vec3(-0.495816,	-0.469658,	-0.730471),
                                         vec3(-0.810407,	-0.203773,	-0.549288),
                                         vec3(0.702667,	0.655321,	0.277152),
                                         vec3(-0.263226,	-0.459216,	0.84843),
                                         vec3(-0.586315,	0.573607,	-0.572022),
                                         vec3(0.597548,	-0.799493,	-0.0612211),
                                         vec3(-0.503628,	-0.694559,	0.513758),
                                         vec3(-0.151531,	-0.9707,	-0.186495),
                                         vec3(-0.406277,	0.670409,	0.620879),
                                         vec3(-0.951893,	0.29543,	0.0813739),
                                         vec3(0.211674,	0.51445,	-0.830985),
                                         vec3(-0.170296,	0.718735,	-0.674106),
                                         vec3(0.934453,	-0.171308,	0.312173));

const vec3 lightColors [64] = vec3[](
                                     vec3(0.00521182,	0.00212318,	0.00237503),
                                     vec3(0.775783,	0.716417,	1.40207),
                                     vec3(0.0622671,	0.0287172,	0.0249731),
                                     vec3(0.012563,	0.00314269,	0.00169478),
                                     vec3(0.00455539,	0.00144987,	0.0012757),
                                     vec3(0.0295916,	0.00649812,	0.00352682),
                                     vec3(0.0249239,	0.00679184,	0.00291705),
                                     vec3(0.00238518,	0.00141129,	0.00204865),
                                     vec3(0.0437176,	0.0144649,	0.017611),
                                     vec3(0.340916,	0.256254,	0.261499),
                                     vec3(0.211026,	0.182618,	0.557771),
                                     vec3(0.0308876,	0.0126275,	0.00894882),
                                     vec3(0.200207,	0.15301,	0.180585),
                                     vec3(0.0201448,	0.00661621,	0.0079783),
                                     vec3(0.0327392,	0.0152856,	0.0110586),
                                     vec3(0.0576628,	0.0309167,	0.0143553),
                                     vec3(0.00734177,	0.00334602,	0.00261892),
                                     vec3(0.0937885,	0.0277603,	0.0208181),
                                     vec3(0.403944,	0.154047,	0.0265064),
                                     vec3(0.0473079,	0.0168941,	0.016979),
                                     vec3(0.0421994,	0.0427472,	0.0583984),
                                     vec3(1.30981,	0.45835,	0.0656782),
                                     vec3(0.121076,	0.0989,	0.151274),
                                     vec3(0.0175308,	0.00488294,	0.00245932),
                                     vec3(0.393723,	0.269764,	0.0719964),
                                     vec3(0.00271467,	0.0014498,	0.00185575),
                                     vec3(0.00488874,	0.00195944,	0.00230665),
                                     vec3(0.010277,	0.00231569,	0.00207359),
                                     vec3(0.0124398,	0.00313595,	0.00305611),
                                     vec3(0.0770417,	0.0177562,	0.00756213),
                                     vec3(0.238802,	0.0353344,	0.00770716),
                                     vec3(0.0170882,	0.00362913,	0.00173228),
                                     vec3(0.0135536,	0.00706163,	0.00871019),
                                     vec3(1.08106,	0.20253,	0.0225603),
                                     vec3(0.488199,	0.38528,	0.511543),
                                     vec3(0.176603,	0.0744748,	0.0663485),
                                     vec3(0.0187279,	0.00515231,	0.00386573),
                                     vec3(0.00436664,	0.00211499,	0.00249297),
                                     vec3(0.0287609,	0.0103246,	0.0127301),
                                     vec3(0.0178167,	0.00300949,	0.00295739),
                                     vec3(0.0190089,	0.00719246,	0.00997383),
                                     vec3(0.00852195,	0.00296741,	0.00310974),
                                     vec3(0.00242544,	0.00137482,	0.00190535),
                                     vec3(0.00304572,	0.00186093,	0.00298483),
                                     vec3(0.00548684,	0.00266476,	0.00320004),
                                     vec3(0.00408536,	0.00182017,	0.00237012),
                                     vec3(0.00568377,	0.00186819,	0.00188315),
                                     vec3(0.00367147,	0.00231524,	0.00395412),
                                     vec3(0.420982,	0.255975,	0.230577),
                                     vec3(0.00643328,	0.00203989,	0.00199598),
                                     vec3(0.0393984,	0.00782152,	0.00500049),
                                     vec3(0.0284206,	0.00976346,	0.0103182),
                                     vec3(0.02676,	0.00813909,	0.00355215),
                                     vec3(0.17268,	0.117015,	0.117033),
                                     vec3(0.153013,	0.0458182,	0.0128394),
                                     vec3(0.0067157,	0.00312631,	0.00347936),
                                     vec3(0.00320873,	0.00182941,	0.0024051),
                                     vec3(0.188296,	0.0560672,	0.0260636),
                                     vec3(0.00278784,	0.00177925,	0.00274015),
                                     vec3(0.0106413,	0.00488356,	0.00582031),
                                     vec3(0.417358,	0.329016,	0.498265),
                                     vec3(0.322286,	0.232472,	0.435887),
                                     vec3(0.639879,	0.601229,	1.19093),
                                     vec3(0.00514545,	0.00200249,	0.0021059));
	
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
// Aufgabe 3.1 

vec3 sphericalToCartesian(vec3 spherical) {
	
  float phi = spherical.x;
  float theta = spherical.y;
  float r = spherical.z;

  vec3 cart;

  cart.x = sin(theta) * cos(phi) * r;
  cart.y = sin(theta) * sin(phi) * r;
  cart.z = cos(theta) * r;
	
  return cart;
}


vec3 cartesianToSpherical(vec3 cart) {
	
  float phi = 0.0;
  float theta = 0.0;
  float r = 0.0;

  r = sqrt(pow(cart.x,2)+pow(cart.y,2)+pow(cart.z,2));
  if (phi != 0){
    phi = atan(cart.y/cart.x);
  }
  else{
    phi = (carty > 0 ? // pi : -pi
  }
  theta = acos(cart.z/r);
	
  return vec3(phi, theta, r);
}

// Ende Aufgabe 3.1 

vec3 sampleEnvMap(vec3 dir)
{

  // Aufgabe 3.3
  vec3 color = vec3( 0.5);
  vec3 angles;
  vec2 coords = vec2(0.5);
	
  color = texture2D( envMapTexture, coords).xyz;

  // Ende Aufgabe 3.3

  return color;
}

vec3 getReflectedDirection(vec3 dir, vec3 normal) 
{
  vec3 reflected = dir - normal*(2.0*dot(normal,dir));

  return reflected;
}

vec3 compProd(vec3 a, vec3 b) {
  return vec3(a.r * b.r, a.g * b.g, a.b * b.b);
}

// n = normal; v = ray.direction; l = light.direction
vec3 BlinnPhongBRDF( vec3 n, vec3 v, vec3 l, vec3 diffuseCol, vec3 specularCol, float sigma) {

  vec3 result = vec3(0.0);
	
  return result;
}

// n = normal; v = ray.direction; l = light.direction
vec3 wardBRDF( vec3 n, vec3 v, vec3 l, vec3 diffuseCol, vec3 specularCol, float sigma) {
  vec3 result = vec3(0.0);
  
  return result;
}



vec4 radiance(Ray ray, vec3 eye)
{
  vec4 color = vec4(0.0, 0.0, 0.0, 1.0);

  float minDistance = 1000000000.0;

  float planeDistance = 0.0;
	
  // ray-g_plane-intersection
  if (false && intersectPlane(ray, g_plane, planeDistance))
	{
      minDistance = planeDistance;
      vec3 p = ray.origin + minDistance * ray.direction;

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
  // no ray-g_plane-intersection --> set background color
  else
	{
      color = vec4(0.2, 0.0, 0.2, 1.0);
      color = vec4( sampleEnvMap( ray.direction), 1.0);
			
	}

  // ray-sphere-intersection
  for (int i = 0; i < NUM_SPHERES; i++)
	{
      float t1, t2;
		
      if (intersectSphere(ray, g_spheres[i], t1, t2) && sphere_switch > 0)
		{
			
          float distance = min(t1, t2);
			
          if (distance > 0.0 && distance < minDistance)
			{
              minDistance = distance;
				
              vec3 p = ray.origin + minDistance * ray.direction;
		
              vec3 n = normalize(p - g_spheres[i].center);

              vec3 l = normalize(vec3(1.0, 5.0, 0.0) - p);
              vec3 v = normalize(eye - p);
              vec3 h = normalize(v + l);

              float angle = dot(n, l);
              float spec = pow(max(0.0, dot(h, n)), 50.0);

              // Aufgabe 3.4
              // reflexion
              vec3 r;
              vec4 refcolor; 

              switch(BRDF_number){
              case(0):
                
				// Ende Aufgabe 3.4
                break;
              case(1):
                // Aufgabe 3.6 Lambert
                color = vec4(0.0, 0.0, 0.0, 1.0);
                for(int j = 0; j < numLights; ++j) 
                  {
                  }
					
                color /= float( numLights);
                // Ende Aufgabe 3.6
                break;
              case(2):
				// Aufgabe 3.7 Ward
                
                color = vec4(0.0, 0.0, 0.0, 1.0);
                for(int j = 0; j < numLights; ++j) 
                  {
				  
                  }
						
                color /= float( numLights);
				// Ende Aufgabe 3.7
                break;
              case(3):
			    // Aufgabe 3.8 Blinn-Phong
                
                color = vec4(0.0, 0.0, 0.0, 1.0);
                for(int j = 0; j < numLights; ++j) 
                  {
                  }
			    // Ende Aufgabe 3.8
                
                color /= float( numLights);
                break;
              case(4):
              default:
                color = g_spheres[i].color;
              }

			}
		}
	}

  color.a = 1.0;
  return color;
}


vec4 gamma(vec4 color)
{
  // Gamma-Korrektur
		
  vec4 amped_col = color;
  if(gain_switch > 0){
    // Aufgabe 3.5
    amped_col.r = color.r*gain_number;
    amped_col.g = color.g*gain_number;
    amped_col.b = color.b*gain_number;
   
    // Ende Aufgabe 3.5
  }
		
  clamp(amped_col.rgb, 0.0, 1.0);

  const float magic = 0.0031308;
  const float factor = 12.92;
  const float a = 0.055;
  const float aPlus1 = a + 1.0;
  const float inv24 = 1.0 / 2.4;

  color.r = (amped_col.r <= magic) ? (amped_col.r * factor) : (aPlus1 * pow(amped_col.r, inv24) - a);
  color.g = (amped_col.g <= magic) ? (amped_col.g * factor) : (aPlus1 * pow(amped_col.g, inv24) - a);
  color.b = (amped_col.b <= magic) ? (amped_col.b * factor) : (aPlus1 * pow(amped_col.b, inv24) - a);
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
  if(BRDF_number == 1){
    for (int i = 0; i < NUM_SPHERES; i++){
      g_spheres[i].color = vec4(1.0, 1.0, 1.0, 1.0);
    }
  }
								
}


void initCamera()
{
  vec3 vert = vec3(0.0, 1.0, 0.0);

  // Aufgabe 3.2
  g_lookat = vec3(0.0, 0.0, 0.0);
  // Ende Aufgabe 3.2
	
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
  return makeRay(g_eye + (((2.0 * x / resolution.x) - 1.0) * g_right) + (((2.0 * y / resolution.y) - 1.0) * g_up)
                 , g_forward
                 );
}


Ray getRayPerspective(float x, float y)
{
	
  return makeRay(g_eye
                 , g_forward + (((2.0 * x / resolution.x) - 1.0) * g_right) + (((2.0 * y / resolution.y) - 1.0) * g_up)
                 );
}


vec4 getColor(float x, float y)
{
  Ray ray = getRayPerspective(x, y);			// Rufen Sie hier getRayPerspective(x, y) auf
	
  vec4 color = radiance(ray, g_eye);

  return color;
}

void main(void)
{
  initWorld();
  initCamera();

  vec4 color = getColor(gl_FragCoord.x, gl_FragCoord.y);
  //color = vec4(0.3,0.5, 1.0 ,0.5);
  color = gamma(color);
	
  frag_color = color;
}
