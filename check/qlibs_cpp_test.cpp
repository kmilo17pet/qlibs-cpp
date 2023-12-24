#include <iostream>
#include <qlibs.h>

using namespace std;

void test_fis( void );
void test_fis2( void );
void test_tdl( void );
void test_fp16( void );
void test_crc( void );

void test_tdl( void )
{
    cout << "TDL test"<< endl; 

    real_t delays[ 10 ];
    tdl delayLine( delays );

    delayLine( 3.5 );
    delayLine( 5 );
    delayLine( 5.0 );
    delayLine( 0.5 );
    delayLine( 1.5 );
    delayLine( 2.5 );
    delayLine( 3.5 );
    delayLine( 4.5 );
    delayLine( 5.5 );
    delayLine( 6.5 );
    delayLine( 7.5 );
    delayLine( 8.5 );

    cout << delayLine[ 0 ] << endl;
    cout << delayLine[ 1 ] << endl;
    cout << delayLine[ -1 ] << endl;
    cout << delayLine[ 100 ] << endl;
}

void test_fis2( void )
{
    enum : fisTag { wt, dax, day, ae };
    /*  I/O Membership functions tags */
    enum : fisTag { phit, thetat };
    enum : fisTag { wtSLOW, wtMED, wtFAST, daxLOW, daxMED, daxHIGH, dayLOW, dayMED, dayHIGH, aeLOW, aeMED, aeHIGH };
    enum : fisTag { phitGYRO, phitBOTH, phitACCEL, thetatGYRO, thetatBOTH, thetatACCEL };

    fis flexnav;
    fisInput flexnav_inputs[4];
    fisOutput flexnav_outputs[2];
    /* I/O Membership Objects */
    fisMF MFin[12], MFout[6];
    static const fisRules rules[] = {
        FIS_RULES_BEGIN
            IF wt IS_NOT wtSLOW THEN phit IS phitGYRO AND thetat IS thetatGYRO END
            IF dax IS daxHIGH THEN thetat IS thetatGYRO END
            IF day IS dayHIGH THEN thetat IS thetatGYRO END
            IF ae IS aeHIGH THEN phit IS phitGYRO AND thetat IS thetatGYRO END
            IF wt IS wtSLOW AND dax IS daxLOW AND ae IS aeLOW THEN phit IS phitACCEL END
            IF wt IS wtSLOW AND day IS dayLOW AND ae IS aeLOW THEN thetat IS thetatACCEL END
            IF wt IS wtSLOW AND dax IS daxLOW AND ae IS aeMED THEN phit IS phitBOTH END
            IF wt IS wtSLOW AND day IS dayLOW AND ae IS aeMED THEN thetat IS thetatBOTH END
            IF wt IS wtSLOW AND dax IS daxMED AND ae IS aeLOW THEN phit IS phitBOTH END
            IF wt IS wtSLOW AND day IS dayMED AND ae IS aeLOW THEN thetat IS thetatBOTH END
            IF wt IS wtMED AND dax IS daxLOW AND ae IS aeLOW THEN phit IS phitBOTH END
            IF wt IS wtMED AND day IS dayLOW AND ae IS aeLOW THEN thetat IS thetatBOTH END
            IF wt IS wtMED AND dax IS_NOT daxLOW THEN phit IS phitGYRO END
            IF wt IS wtMED AND day IS_NOT dayLOW THEN thetat IS thetatGYRO END
            IF wt IS wtMED AND ae IS_NOT aeLOW THEN phit IS phitGYRO AND thetat IS thetatGYRO END
        FIS_RULES_END
    };
    real_t rulesStrength[ 15 ] = { 0.0 };
    //flexnav.setup( specs );
    flexnav.setup( Mamdani, flexnav_inputs, flexnav_outputs, MFin, MFout, rules, rulesStrength );
    flexnav.setupInput( wt, 0.0, 0.5 );
    flexnav.setupInput( dax, 0.0, 5.0 );
    flexnav.setupInput( day, 0.0, 5.0 );
    flexnav.setupInput( ae, 0.0, 20.0 );
    flexnav.setupOutput( phit, 0.0, 1.0 );
    flexnav.setupOutput( thetat, 0.0, 1.0 );
    flexnav.setInputMF( wt, wtSLOW, trimf, (const real_t []){ -0.2 ,0.0 ,0.2 } );
    flexnav.setInputMF( wt, wtMED, trimf, (const real_t []){ 0.1 ,0.25 , 0.4 } );
    flexnav.setInputMF( wt, wtFAST, trimf, (const real_t []){ 0.3 ,0.5 ,0.7 });
    flexnav.setInputMF( dax, daxLOW, trimf, (const real_t []){ -1.0 ,0.0 ,2.0 } );
    flexnav.setInputMF( dax, daxMED, trimf, (const real_t []){ 1.0 , 2.5 , 4.0 } );
    flexnav.setInputMF( dax, daxHIGH, trimf, (const real_t []){ 3.0 ,5.0 , 7.0 } );
    flexnav.setInputMF( day, dayLOW, trimf, (const real_t []){ -2.0 ,0.0 , 2.0 } );
    flexnav.setInputMF( day, dayMED, trimf, (const real_t []){ 1.0 ,2.5 ,4.0 } );
    flexnav.setInputMF( day, dayHIGH, trimf, (const real_t []){ 3.0 ,5.0 , 7.0 } );
    flexnav.setInputMF( ae, aeLOW, trimf, (const real_t []){ -8.0 ,0.0 ,8.0 } );
    flexnav.setInputMF( ae, aeMED, trimf, (const real_t []){ 5.0 ,10.0 , 15.0 } );
    flexnav.setInputMF( ae, aeHIGH, trimf, (const real_t []){ 12.0 ,20.0 ,28.0 } );
    
    flexnav.setOutputMF( phit, phitGYRO, trimf, (const real_t []){ -0.4 ,0.0 ,0.4 } );
    flexnav.setOutputMF( phit, phitBOTH, trimf, (const real_t []){ 0.2 , 0.5 , 0.8 } );
    flexnav.setOutputMF( phit, phitACCEL, trimf, (const real_t []){ 0.6 , 1.0 , 1.4 } );
    flexnav.setOutputMF( thetat, thetatGYRO, trimf, (const real_t []){ -0.4 ,0.0 , 0.4 } );
    flexnav.setOutputMF( thetat, thetatBOTH, trimf, (const real_t []){ 0.2 , 0.5 , 0.8 } );
    flexnav.setOutputMF( thetat, thetatACCEL, trimf, (const real_t []){ 0.6 , 1.0 , 1.4 } );

    flexnav << wt << 0 << dax << 0 << day << 3 << ae << 0;
    flexnav.fuzzify();
    if ( flexnav.inference() ) {
        flexnav.deFuzzify();
    }
    cout << "phit = " << flexnav[ phit ] << " thetat = " << flexnav[ thetat ] << endl;

}

void test_fis( void )
{
    real_t xag[ 100 ], yag[ 100 ];
    // I/O Names
    enum : fisTag { service, food};
    enum : fisTag { tip};
    // I/O Membership functions tags
    enum : fisTag { poor, good, excellent, rancid, delicious };
    enum : fisTag { cheap, average, generous};

    fis tipper;
    // I/O Fuzzy Objects
    fisInput tipper_inputs[ 2 ];
    fisOutput tipper_outputs[ 1 ];
    // I/O Membership Objects
    fisMF MFin[5], MFout[3];

    const fisRules rules[] = { 
        FIS_RULES_BEGIN
            IF service IS poor OR food IS rancid THEN tip IS cheap END
            IF service IS good THEN tip IS average END
            IF service IS excellent OR food IS delicious THEN tip IS generous END
        FIS_RULES_END
    };
    real_t rulesStrength[ 3 ];

    tipper.setup( Mamdani, tipper_inputs, tipper_outputs, MFin, MFout, rules, rulesStrength );
    tipper.setupInput( service, 0.0, 1.0 );
    tipper.setupInput( food, 0.0, 10.0 );
    tipper.setupOutput( tip, 0.0, 30.0 );
    tipper.setInputMF( service, poor, gaussmf, (const real_t[]){ 1.5, 0.0 } );
    tipper.setInputMF( service, good, gaussmf, (const real_t[]){ 1.5, 5.0 } );
    tipper.setInputMF( service, excellent, gaussmf, (const real_t[]){ 1.5, 10.0 } );
    tipper.setInputMF( food, rancid, trapmf, (const real_t[]){ 0.0, 0.0, 1.0, 3.0 } );
    tipper.setInputMF( food, delicious, trapmf, (const real_t[]){ 7.0, 9.0, 10.0, 10.0 } );
    tipper.setOutputMF( tip, cheap, trimf, (const real_t[]){ 0.0, 5.0, 10.0 } );
    tipper.setOutputMF( tip, average, trimf, (const real_t[]){10.0, 15.0, 20.0 } );
    tipper.setOutputMF( tip, generous, trimf, (const real_t[]){ 20.0, 25.0, 30.0 } );
    tipper_outputs[ tip ].storeAggregatedRegion( xag, yag );


    tipper << service << 1 << food << 9;
    tipper.fuzzify();
    if ( tipper.inference() ) {
        tipper.deFuzzify();
    }
    cout << "tip_out = " << tipper[ tip ] << endl;
}

void test_fp16( void ) 
{
    fp16 x = 6.5_fp;
    fp16 y = 7.33_fp;
    fp16 z = x * y;
    fp16 w = z/2.0_fp;
    fp16 r = -y;
    fp16 t = fp16::exp( 0.1_fp );
    cout << x.raw() << " "<< y.raw() << " " << z.raw() << endl;
    cout <<  z  << endl;
    cout <<  fp16::sqrt( z )  << endl;
    cout <<  w  << endl;
    cout <<  r  << endl;
    cout <<  t  << endl;
    cout <<  2.0_fp * ( 5.3_fp - 2.1_fp )  << endl;
    cout <<  2_fp / 3.5_fp  << endl;
    cout <<  fp16::log( 3.5_fp )  << endl;
    cout <<  ( 2_fp*( x + y ) )/y  << endl;
    cout <<  x  << endl;
    x += 1.3_fp;
    cout <<  x  << endl;
    x -= 0.5_fp;
    cout <<  x << endl;
    x += y;
    cout <<  x  << endl;
    x++;
    cout << x << endl;
}

void test_crc( void )
{
    crc crcx;
    cout << sizeof(crc) << endl;
}

int main()
{
    test_crc();
    test_fp16();
    return EXIT_SUCCESS;
    cout << "continuousSystem"<< endl; 
    constexpr size_t SYS_ORDER = 3;
    real_t num[] = { 0.0, 2.0, 3.0, 6.0 };
    real_t den[] = { 1.0, 6.0, 11.0, 16.0 };
    continuousStates xC = { 0, 0, 0 }; 
    continuousSystem gc( num, den, xC, SYS_ORDER+1, 0.01 );

    for ( int i = 0; i < 3005; i++ ) {
        cout << gc.excite( 1.0 ) << endl;
    }

    cout << "discreteSystem"<< endl; 
    constexpr size_t NB = 3;
    constexpr size_t NA = 3;
    real_t nd[] = { 0.1, 0.2, 0.3 };      // 0.1 + 0.2*z^(-1) + 0.3*z^(-2)
    real_t dd[] = { 1.0, -0.85, 0.02 };   // 1 - 0.85*z^(-1) + 0.02*z^(-2)
    discreteStates xD = { 0, 0, 0};
    discreteSystem gd( nd ,dd, xD, NB, NA );
    for ( int i = 0; i < 35; i++ ) {
        cout << gd.excite( 1.0 ) << endl;
    }



    test_tdl();
    test_fis();
    test_fis2();
    return 0;
}