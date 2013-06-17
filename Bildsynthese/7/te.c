#version 420

uniform mat4 modelviewprojection_matrix;
uniform mat3 normal_matrix;
uniform int time;

layout(quads, equal_spacing, cw) in;
out vec3 tess_normal;

#define M_PI 3.1415926535897932384626433832795

void main() {
    float u = gl_TessCoord.x;
    float v = gl_TessCoord.y;

    // Bilineare Interpolation
    vec3 p0 = mix(gl_in[0].gl_Position.xyz, gl_in[1].gl_Position.xyz, u);
    vec3 p1 = mix(gl_in[3].gl_Position.xyz, gl_in[2].gl_Position.xyz, u);
    vec3 p2 = mix(p0, p1, v);
    
    // Displacement
    p2.z = sin(u*M_PI) * sin(v*M_PI);
    
    // Animation
    float step = mod(time,1000);
    p2.z *= cos(step);
    
    gl_Position = modelviewprojection_matrix * vec4(p2, 1.0);

    // Normal 
    vec3 grad = vec3(1.0);
    grad.x = cos(u*M_PI) * M_PI * sin(v*M_PI);
    grad.y = sin(u*M_PI) * cos(v*M_PI) * M_PI;
    grad.z = cos(step);
    tess_normal = grad;
}

