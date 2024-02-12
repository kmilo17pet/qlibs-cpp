#include <qlibs.h>

smootherGMWF filter;
constexpr real_t SMOOTHER_SAMPLE_TIME = 20; /*mS*/
int inputPin = A0;   // select the input pin for the potentiometer that will drive the system

void setup() {
  Serial.begin(9600);
  /*setup the gaussian filter*/
  constexpr size_t SMOOTHER_WINDOW_SIZE = 30;
  constexpr real_t GAUSS_SIGMA = 20.0f;
  constexpr real_t GAUSS_CENTER = 0.0f;
  real_t window[ SMOOTHER_WINDOW_SIZE ] = { 0.0f };
  real_t kernel[ SMOOTHER_WINDOW_SIZE ] = { 0.0f };
  filter.setup( GAUSS_SIGMA, GAUSS_CENTER, window, kernel );
}

void loop() {
  real_t noise = static_cast<float>( random(-1000, 1000) )/10000.00;
  real_t noised_signal =  noise + map( analogRead( inputPin ), 0, 1023, 0.0f, 100.0f );
  real_t filterd_signal = filter.smooth( noised_signal );
  Serial.print( noised_signal );
  Serial.print( "," );
  Serial.println( filterd_signal );
  delay( SMOOTHER_SAMPLE_TIME);
}
