#version 330

uniform vec2 resolution;
uniform float timer;

//output color
layout(location = 0) out vec4 gl_FragColor;

//incoming coordinates in range [0 ..1]
in vec2 texCoords;

// here you can compute the the color for a point with texture coordinates texCoords
// the output color must be in gl_FragColor : gl_FragColor = vec4(R, G, B, 1.0f);
// with resolution.x you can access the current window width
// with resolution.y you can access the current window height

//with texCoords.x you can access the x component
//with texCoords.y you can access the y component

void main( void )
{
	//insert your code here
    vec4 blue = vec4(0.0f, 0.0f, 1.0f, 1.0f);
    vec4 white = vec4(1.0f, 1.0f, 1.0f, 1.0f);
    vec4 red = vec4(1.0f, 0.0f, 0.0f, 1.0f);
    vec4 color = vec4(1.0f, 0.0f, 0.0f, 1.0f);
    
    float alpha = timer;
    float beta = resolution.x;
    float gamma = resolution.y;
    float phi1 = texCoords.x;
    float phi2 = texCoords.y;
    
	float value = abs((sin(alpha)+1.0)*sin(beta*phi1+cos(alpha))*cos(gamma*phi2+sin(alpha)));
    
    if ( value >= 0.0f && value <= 0.5f)
        color = mix(blue, white, (value / 0.5f));
    else if( value > 0.5f && value <= 1.0f)
        color = mix(white, red, ((value-0.5f)/0.5f));
	
	gl_FragColor = color;
}

