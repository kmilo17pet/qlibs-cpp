#include <qlibs.h>

using namespace qlibs;

continuousTF<3> ctf= {
    { 0.0f, 2.0f, 3.0f, 6.0f },
    { 1.0f, 6.0f, 11.0f, 16.0f },
};
real_t dt = 0.05f;
continuousSystem gc( ctf, dt );

int inputPin = A0;   // select the input pin for the potentiometer that will drive the system

void setup() {
  Serial.begin(9600);
}

void loop() {
  real_t ut = map( analogRead( inputPin ), 0, 1023, 0.0f, 100.0f );
  real_t yt = gc.excite( ut );
  Serial.print( ut );
  Serial.print( "," );
  Serial.println( yt );
  delay( dt*1000 );
}
