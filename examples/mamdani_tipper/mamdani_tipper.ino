#include <qlibs.h>

enum : fis::tag { service, food};
enum : fis::tag { tip};
// I/O Membership functions tags
enum : fis::tag { poor, good, excellent, rancid, delicious };
enum : fis::tag { cheap, average, generous};
fis::instance tipper;
// I/O Fuzzy Objects
fis::input tipper_inputs[ 2 ];
fis::output tipper_outputs[ 1 ];
// I/O Membership Objects
fis::mf MFin[5], MFout[3];
const fis::rules rules[] = {
    FIS_RULES_BEGIN
        IF service IS poor OR food IS rancid THEN tip IS cheap END
        IF service IS good THEN tip IS average END
        IF service IS excellent OR food IS delicious THEN tip IS generous END
    FIS_RULES_END
};
real_t rulesStrength[ 3 ];

int foodPin = A0;   // select the input pin for the potentiometer that drives the food score
int servicePin = A1;   // select the input pin for the potentiometer that drives the service score

void setup() {
  Serial.begin( 9600 );
  tipper.setup( fis::Mamdani, tipper_inputs, tipper_outputs, MFin, MFout, rules, rulesStrength );
  tipper.setupInput( service, 0.0f, 1.0f );
  tipper.setupInput( food, 0.0f, 10.0f );
  tipper.setupOutput( tip, 0.0f, 30.0f );
  tipper.setupInputMF( service, poor, fis::gaussmf, (const real_t[]){ 1.5f, 0.0f } );
  tipper.setupInputMF( service, good, fis::gaussmf, (const real_t[]){ 1.5f, 5.0f } );
  tipper.setupInputMF( service, excellent, fis::gaussmf, (const real_t[]){ 1.5f, 10.0f } );
  tipper.setupInputMF( food, rancid, fis::trapmf, (const real_t[]){ 0.0f, 0.0f, 1.0f, 3.0f } );
  tipper.setupInputMF( food, delicious, fis::trapmf, (const real_t[]){ 7.0f, 9.0, 10.0f, 10.0f } );
  tipper.setupOutputMF( tip, cheap, fis::trimf, (const real_t[]){ 0.0f, 5.0f, 10.0f } );
  tipper.setupOutputMF( tip, average, fis::trimf, (const real_t[]){10.0f, 15.0f, 20.0f } );
  tipper.setupOutputMF( tip, generous, fis::trimf, (const real_t[]){ 20.0f, 25.0f, 30.0f } );
}

void loop() {
  real_t serviceScoreValue = map( analogRead(servicePin), 0, 1023, 0.0f, 1.0f );
  real_t foodScoreValue = map( analogRead(foodPin), 0, 1023, 0.0f, 10.0f );

  tipper.setInput( service, serviceScoreValue );
  tipper.setInput( food, foodScoreValue );

  tipper.fuzzify();
  if ( tipper.inference() ) {
      tipper.deFuzzify();
  }
  Serial.print("food[in] = ");
  Serial.println(foodScoreValue);
  Serial.print("service[in] = ");
  Serial.println(serviceScoreValue);
  Serial.print("tip[out] = ");
  Serial.println( tipper[ tip ] );
  delay( 250 );
}
