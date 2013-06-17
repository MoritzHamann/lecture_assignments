#version 420

uniform mat4 modelviewprojection_matrix;

layout(location = 0) in vec3 in_position;
out vec3 vPosition;


void main() {
    gl_Position = vec4(in_position.xyz, 1.0);
    vPosition.xyz = in_position.xyz;
}
