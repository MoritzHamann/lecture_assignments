#version 420

#define M_PI 3.1415926535897932384626433832795

layout(location = 0) in vec3 in_position;
uniform mat4  modelviewprojection_matrix;
uniform int time;
uniform int speed;
//out vec3 vPosition;



// TODO Aufgabe 8.2.1
float timeperiod(int time, int speed){
  // divide speed to get slower movement
  float offset = 0.0;
  offset = 0.5 * sin(speed/10.00 * time/1000.00) + 1;
  return offset;
}
// ENDE Aufgabe 8.2.1


void main() {

  // TODO Aufgabe 8.2.2
  vec3 pos = in_position.xyz;
  float itime = timeperiod(time, speed);
  pos.x += itime;
  pos.y += itime;
  pos.z += itime;
  // ENDE Aufgabe 8.2.2

  gl_Position = modelviewprojection_matrix*vec4(pos, 1.0);
}
