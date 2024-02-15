#include <qlibs.h>

continuousTF<1> processTransferFunction= {
    { 0.0f, 1.5f, },
    { 2.0f, 1.0f, },
};

constexpr real_t dt = 0.05f;    /*Time step*/
real_t yt = 0.0f;               /*Process output*/
real_t ut = 0.0f;               /*Controller output*/
real_t wt = 0.0f;               /*Set-Point*/
continuousSystem process( processTransferFunction, dt );
pidController controller;
int inputPin = A0;   // select the input pin for the potentiometer that will drive the setPoint

void setup() {
  Serial.begin(9600);
  controller.setup( 5.0f, 10.0f ,0.01f ,dt );
  controller.setSaturation( 0.0f, 100.0f );
}

void loop() {
  /*Set-Point is driven from a potentiometer connected to analog input*/
  wt = map( analogRead( inputPin ), 0, 1023, 0.0f, 100.0f );
  ut = controller.control( wt, yt );
  yt = process.excite( ut );

  Serial.print( ut );
  Serial.print( "," );
  Serial.print( wt );
  Serial.print( "," );
  Serial.println( yt );
  delay( dt*1000 );
}
