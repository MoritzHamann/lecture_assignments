#version 330

layout(location = 0) in vec4  in_position;

uniform mat4 mvp;

out vec2 texCoords;

void main()
{
    vec4 vert = in_position; 
    gl_Position = mvp * vert; 
    texCoords = in_position.xy;
}