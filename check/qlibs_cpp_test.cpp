#include <iostream>
#include <qlibs.h>

using namespace std;

void test_fis( void );
void test_fis2( void );
void test_tdl( void );
void test_fp16( void );
void test_crc( void );
void test_ltisys( void );

void test_tdl( void )
{
    cout << "TDL TEST"<< endl; 

    real_t delays[ 10 ];
    tdl delayLine( delays );

    delayLine( 3.5f );
    delayLine( 5.0f );
    delayLine( 5.0f );
    delayLine( 0.5f );
    delayLine( 1.5f );
    delayLine( 2.5f );
    delayLine( 3.5f );
    delayLine( 4.5f );
    delayLine( 5.5f );
    delayLine( 6.5f );
    delayLine( 7.5f );
    delayLine( 8.5f );

    cout << delayLine[ 0 ] << endl;
    cout << delayLine[ 1 ] << endl;
    cout << delayLine[ -1 ] << endl;
    cout << delayLine[ 100 ] << endl;
}

void test_fis2( void )
{
    cout << "FIS2 TEST"<< endl; 
    enum : fisTag { wt, dax, day, ae };
    enum : fisTag { phit, thetat };
    enum : fisTag { wtSLOW, wtMED, wtFAST, daxLOW, daxMED, daxHIGH, dayLOW, dayMED, dayHIGH, aeLOW, aeMED, aeHIGH };
    enum : fisTag { phitGYRO, phitBOTH, phitACCEL, thetatGYRO, thetatBOTH, thetatACCEL };

    fisSystem<Mamdani,4, 2, 12, 6, 15> flexnav = {
        (const fisRules[]) {
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
        }
    };

    flexnav.setup();
    flexnav.setupInput( wt, 0.0f, 0.5f );
    flexnav.setupInput( dax, 0.0f, 5.0f );
    flexnav.setupInput( day, 0.0f, 5.0f );
    flexnav.setupInput( ae, 0.0f, 20.0f );
    flexnav.setupOutput( phit, 0.0f, 1.0f );
    flexnav.setupOutput( thetat, 0.0f, 1.0f );
    flexnav.setInputMF( wt, wtSLOW, trimf, (const real_t []){ -0.2f ,0.0f ,0.2f } );
    flexnav.setInputMF( wt, wtMED, trimf, (const real_t []){ 0.1f ,0.25f , 0.4f } );
    flexnav.setInputMF( wt, wtFAST, trimf, (const real_t []){ 0.3f ,0.5f ,0.7f });
    flexnav.setInputMF( dax, daxLOW, trimf, (const real_t []){ -1.0f ,0.0f ,2.0f } );
    flexnav.setInputMF( dax, daxMED, trimf, (const real_t []){ 1.0f , 2.5f , 4.0f } );
    flexnav.setInputMF( dax, daxHIGH, trimf, (const real_t []){ 3.0f ,5.0f , 7.0f } );
    flexnav.setInputMF( day, dayLOW, trimf, (const real_t []){ -2.0f ,0.0f , 2.0f } );
    flexnav.setInputMF( day, dayMED, trimf, (const real_t []){ 1.0f ,2.5f ,4.0f } );
    flexnav.setInputMF( day, dayHIGH, trimf, (const real_t []){ 3.0f ,5.0f , 7.0f } );
    flexnav.setInputMF( ae, aeLOW, trimf, (const real_t []){ -8.0f ,0.0f ,8.0f } );
    flexnav.setInputMF( ae, aeMED, trimf, (const real_t []){ 5.0f ,10.0f , 15.0f } );
    flexnav.setInputMF( ae, aeHIGH, trimf, (const real_t []){ 12.0f ,20.0f ,28.0f } );
    
    flexnav.setOutputMF( phit, phitGYRO, trimf, (const real_t []){ -0.4f ,0.0f ,0.4f } );
    flexnav.setOutputMF( phit, phitBOTH, trimf, (const real_t []){ 0.2f , 0.5f , 0.8f } );
    flexnav.setOutputMF( phit, phitACCEL, trimf, (const real_t []){ 0.6f , 1.0f , 1.4f } );
    flexnav.setOutputMF( thetat, thetatGYRO, trimf, (const real_t []){ -0.4f ,0.0f , 0.4f } );
    flexnav.setOutputMF( thetat, thetatBOTH, trimf, (const real_t []){ 0.2f , 0.5f, 0.8f } );
    flexnav.setOutputMF( thetat, thetatACCEL, trimf, (const real_t []){ 0.6f , 1.0f , 1.4f } );

    flexnav << wt << 0 << dax << 0 << day << 3 << ae << 0;
    flexnav.fuzzify();
    if ( flexnav.inference() ) {
        flexnav.deFuzzify();
    }
    cout << "phit = " << flexnav[ phit ] << " thetat = " << flexnav[ thetat ] << endl;

}

void test_fis( void )
{
    cout << "FIS TEST"<< endl; 
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
    tipper.setupInput( service, 0.0f, 1.0f );
    tipper.setupInput( food, 0.0f, 10.0f );
    tipper.setupOutput( tip, 0.0f, 30.0f );
    tipper.setInputMF( service, poor, gaussmf, (const real_t[]){ 1.5f, 0.0f } );
    tipper.setInputMF( service, good, gaussmf, (const real_t[]){ 1.5f, 5.0f } );
    tipper.setInputMF( service, excellent, gaussmf, (const real_t[]){ 1.5f, 10.0f } );
    tipper.setInputMF( food, rancid, trapmf, (const real_t[]){ 0.0f, 0.0f, 1.0f, 3.0f } );
    tipper.setInputMF( food, delicious, trapmf, (const real_t[]){ 7.0f, 9.0, 10.0f, 10.0f } );
    tipper.setOutputMF( tip, cheap, trimf, (const real_t[]){ 0.0f, 5.0f, 10.0f } );
    tipper.setOutputMF( tip, average, trimf, (const real_t[]){10.0f, 15.0f, 20.0f } );
    tipper.setOutputMF( tip, generous, trimf, (const real_t[]){ 20.0f, 25.0f, 30.0f } );
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
    cout << "FP16 TEST "<<endl;
    fp16 x = 6.5_fp;
    fp16 y = 7.33_fp;
    fp16 z = x * y;
    fp16 w = z/2.0_fp;
    fp16 r = -y;
    fp16 t = fp16::exp( 0.1_fp );
    cout << sizeof( fp16 ) << endl;
    cout << t << endl;
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
    cout << x << endl;
    x++;
    cout << x << endl;
    x = 20.77;
    cout << x << endl;
    fp16 ab = fp16::from( -57.22 );
    cout << ab << endl;
    cout << FP_PI << endl;

    fp16 a = 1.5_fp;
    fp16 b = 5.2_fp;
    fp16 c = 4.0_fp;
    fp16 result;
    result = ( -b + fp16::sqrt( ( b*b )  - ( 4_fp*a*c ) ) )/( 2_fp*a );
    std::cout << " result = " << result << std::endl;
}

void test_crc( void )
{
    cout << "CRC TEST"<< endl; 
    auto res = crc::crc16_A( "hello world", 11 );
    cout << sizeof(crc) << " res = "<< res << endl;
}

void test_ltisys( void )
{
    cout << "LTISYS TEST"<< endl; 
    cout << "continuousSystem"<< endl; 
    continuousTF<3> ctf= {
        { 0.0f, 2.0f, 3.0f, 6.0f },
        { 1.0f, 6.0f, 11.0f, 16.0f },
    };
    continuousSystem gc( ctf, 0.01f );
    for ( int i = 0; i < 1000; i++ ) {
        cout << gc.excite( 1.0f ) << endl;
    }

    cout << "discreteSystem"<< endl; 
    //discreteTF<3,3> dtf= {
    //    { 0.1f, 0.2f, 0.3f },
    //    { 1.0f, -0.85f, 0.02f },
    //};
    real_t num[] = { 0.1f, 0.2f, 0.3f };
    real_t den[] = { 1.0f, -0.85f, 0.02f };
    discreteStates<3> xd= { 0.0f, 0.0f, 0.0f };
    discreteSystem gd( num, den, xd );
    for ( int i = 0; i < 20; i++ ) {
        cout << gd.excite( 1.0f ) << endl;
    }

}

void test_ffmath(void)
{
    cout << "ffmath" << endl;
    cout << ffmath::exp( 0.5f ) << endl;
    cout << ffmath::log( 0.5f ) << endl;
    cout << ffmath::sin( 0.5f ) << endl;
    cout << ffmath::cos( 0.5f ) << endl;
    cout << ffmath::tan( 0.5f ) << endl;
    cout << ffmath::asin( 0.5f ) << endl;
    cout << ffmath::acos( 0.5f ) << endl;
    cout << ffmath::atan( 0.5f ) << endl;
    cout << ffmath::cosh( 0.5f ) << endl;
    cout << ffmath::sinh( 0.5f ) << endl;
    cout << ffmath::tanh( 0.5f ) << endl;
    cout << ffmath::atan2( 0.1f, 0.2f ) << endl;
    cout << ffmath::pow( 3.8f, 2.5f ) << endl;
    cout << ffmath::getNan() << endl;
    cout << ffmath::getInf() << endl;
}

int main()
{
    test_crc();
    test_fp16();
    test_ltisys();
    test_tdl();
    test_fis();
    test_fis2();
    test_ffmath();
    return 0;
}