#version 420
// TODO Aufgabe 8.1 
layout(triangle_strip, max_vertices = 3) out;


void main()
{
    // TODO create primitive
    gl_Position = gl_in[0].gl_Position;
    EmitVertex();
    
    gl_Position = gl_in[1].gl_Position;
    EmitVertex();
    
    gl_Position = gl_in[2].gl_Position;
    EmitVertex();
        
    EndPrimitive();

}  
