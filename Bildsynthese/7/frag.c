#version 420

out vec4 frag_color;
in vec3 tess_normal;

void main() {

	vec3 light = vec3( 0.0, 10.0, 10.0);
	float theta = dot(normalize(tess_normal), normalize(light));
	
    frag_color = vec4(1,0,0,1)* theta;
}

