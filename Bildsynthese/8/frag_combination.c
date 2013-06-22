#version 420

uniform sampler2D timeTexture[30];

uniform int MotionNumber;
uniform int Fade;
uniform int Strong;

layout(location = 0) out vec4 frag_color;

in vec2 texCoords;

void main() {
    vec4 color = vec4(0,0,0,1);
    int i;
    float weights = 0.0;
    vec4 temp_color = vec4(0,0,0,1);
    float factor = 0;
    // TODO Aufgabe 8.3 Motion Blur
    // dies hier ersetzen
    
    for (int i = 0; i < MotionNumber ; i++){
        temp_color = texture2D(timeTexture[i], texCoords);
        if (Fade == 1){
            temp_color.rgb *= 1-(float(i)/MotionNumber);
        }
        
        color += temp_color;
        if(i == 0 && Strong == 1){
            // only make white fragment red
            if (temp_color.r == 1.0 && temp_color.g == 1.0 && temp_color.b == 1.0){
                color = vec4(1.0, 0.0, 0.0, 1.0);
                break; // exit for loop for this fragment since red should not be overdrawn
            }
        }
    }
    frag_color = color;

}
