#version 420

uniform int tessLevelInner;
uniform int tessLevelOuter;
uniform mat4 modelviewprojection_matrix;
uniform mat3 normal_matrix;


layout(vertices = 4) out;
#define ID gl_InvocationID

void main() {

    gl_out[gl_InvocationID].gl_Position = gl_in[gl_InvocationID].gl_Position;
    
    if (gl_InvocationID == 0) {
        gl_TessLevelInner[0] = tessLevelInner;
        gl_TessLevelInner[1] = tessLevelInner;
        gl_TessLevelOuter[0] = tessLevelOuter;
        gl_TessLevelOuter[1] = tessLevelOuter;
        gl_TessLevelOuter[2] = tessLevelOuter;
        gl_TessLevelOuter[3] = tessLevelOuter;
    }

}
