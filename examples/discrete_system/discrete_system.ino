#include <qlibs.h>

using namespace qlibs;

discreteTF<3,3> dtf= {
    { 0.1f, 0.2f, 0.3f },
    { 1.0f, -0.85f, 0.02f },
};
real_t Tm = 0.1f;
discreteSystem gc( dtf );

int inputPin = A0;   // select the input pin for the potentiometer that will drive the system

void setup() {
  Serial.begin(9600);
}

void loop() {
  real_t uk = map( analogRead( inputPin ), 0, 1023, 0.0f, 100.0f );
  real_t yk = gc.excite( uk );
  Serial.print( uk );
  Serial.print( "," );
  Serial.println( yk );
  delay( Tm*1000 );
}
