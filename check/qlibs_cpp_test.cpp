#include <iostream>
#include <qlibs.h>
#include <cmath>
#include <include/mat.hpp>

using namespace std;

void test_fis( void );
void test_fis2( void );
void test_fis3( void );
void test_tdl( void );
void test_fp16( void );
void test_crc( void );
void test_ltisys( void );
void test_ffmath( void );
void test_mat( void );
void test_interp1( void );

void test_fis3( void )
{
     cout << "FIS2 TEST3"<< endl;

    enum { input1, input2 };
    enum { output1, output2 };
    /* I/O Membership functions tags */
    enum { input1_Neg, input1_Pos, input2_Small, input2_Big };
    enum { output1_mf1, output2_mf1, output2_mf2, output2_mf3, output2_mf4 };

    fis::system<fis::Sugeno,2, 2, 4, 5, 4> takaji = {
        (const fis::rules[]) {
            FIS_RULES_BEGIN
                IF input1 IS input1_Neg AND input2 IS input2_Small THEN output1 IS output1_mf1 AND output2 IS output2_mf4 END
                IF input1 IS input1_Neg AND input2 IS input2_Big THEN output1 IS output1_mf1 AND output2 IS output2_mf3 END
                IF input1 IS input1_Pos AND input2 IS input2_Small THEN output1 IS output1_mf1 AND output2 IS output2_mf2 END
                IF input1 IS input1_Pos AND input2 IS input2_Big THEN output1 IS output1_mf1 AND output2 IS output2_mf1 END
            FIS_RULES_END
        }
    };

    takaji.setup();
    takaji.setupInput( input1, 0.5000f, 3.5000f );
    takaji.setupInput( input2, -1.0000f, 4.0000f );
    takaji.setupOutput( output1, 0.0000f, 1.0000f );
    takaji.setupOutput( output2, 0.0000f, 1.0000f );

    /* Add membership functions to the inputs */
    takaji.setupInputMF( input1, input1_Neg, fis::trimf, (const real_t []){ 0.5000f, 0.5000f, 3.5000f } );
    takaji.setupInputMF( input1, input1_Pos, fis::trimf, (const real_t []){ 0.5000f, 3.5000f, 3.5000f } );
    takaji.setupInputMF( input2, input2_Small, fis::trimf, (const real_t []){ -1.0000f, -1.0000f, 4.0000f} );
    takaji.setupInputMF( input2, input2_Big, fis::trimf, (const real_t []){ -1.0000f, 4.0000f, 4.0000f } );
    /* Add membership functions to the outputs */
    takaji.setupOutputMF( output1, output1_mf1, fis::linearmf, (const real_t []){ 0.0000f, 1.0000f, 0.0000f } );
    takaji.setupOutputMF( output2, output2_mf1, fis::linearmf, (const real_t []){ 3.5000f, 4.0000f, 0.0000f } );
    takaji.setupOutputMF( output2, output2_mf2, fis::linearmf, (const real_t []){ 3.5000f, -1.0000f, 0.0000f });
    takaji.setupOutputMF( output2, output2_mf3, fis::linearmf, (const real_t []){ 0.5000f, 4.0000f, 0.0000f } );
    takaji.setupOutputMF( output2, output2_mf4, fis::linearmf, (const real_t []){ 0.5000f, -1.0000f, 0.0000f } );

    takaji << input1 << 0 << input2;
    takaji.fuzzify();
    if ( takaji.inference() ) {
        takaji.deFuzzify();
    }
    cout << "output1 = " << takaji[ output1 ] << " output2 = " << takaji[ output2 ] << endl;
}

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
    enum : fis::tag { wt, dax, day, ae };
    enum : fis::tag { phit, thetat };
    enum : fis::tag { wtSLOW, wtMED, wtFAST, daxLOW, daxMED, daxHIGH, dayLOW, dayMED, dayHIGH, aeLOW, aeMED, aeHIGH };
    enum : fis::tag { phitGYRO, phitBOTH, phitACCEL, thetatGYRO, thetatBOTH, thetatACCEL };

    fis::system<fis::Mamdani,4, 2, 12, 6, 15> flexnav = {
        (const fis::rules[]) {
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
    flexnav.setupInputMF( wt, wtSLOW, fis::trimf, (const real_t []){ -0.2f ,0.0f ,0.2f } );
    flexnav.setupInputMF( wt, wtMED, fis::trimf, (const real_t []){ 0.1f ,0.25f , 0.4f } );
    flexnav.setupInputMF( wt, wtFAST, fis::trimf, (const real_t []){ 0.3f ,0.5f ,0.7f });
    flexnav.setupInputMF( dax, daxLOW, fis::trimf, (const real_t []){ -1.0f ,0.0f ,2.0f } );
    flexnav.setupInputMF( dax, daxMED, fis::trimf, (const real_t []){ 1.0f , 2.5f , 4.0f } );
    flexnav.setupInputMF( dax, daxHIGH, fis::trimf, (const real_t []){ 3.0f ,5.0f , 7.0f } );
    flexnav.setupInputMF( day, dayLOW, fis::trimf, (const real_t []){ -2.0f ,0.0f , 2.0f } );
    flexnav.setupInputMF( day, dayMED, fis::trimf, (const real_t []){ 1.0f ,2.5f ,4.0f } );
    flexnav.setupInputMF( day, dayHIGH, fis::trimf, (const real_t []){ 3.0f ,5.0f , 7.0f } );
    flexnav.setupInputMF( ae, aeLOW, fis::trimf, (const real_t []){ -8.0f ,0.0f ,8.0f } );
    flexnav.setupInputMF( ae, aeMED, fis::trimf, (const real_t []){ 5.0f ,10.0f , 15.0f } );
    flexnav.setupInputMF( ae, aeHIGH, fis::trimf, (const real_t []){ 12.0f ,20.0f ,28.0f } );

    flexnav.setupOutputMF( phit, phitGYRO, fis::trimf, (const real_t []){ -0.4f ,0.0f ,0.4f } );
    flexnav.setupOutputMF( phit, phitBOTH, fis::trimf, (const real_t []){ 0.2f , 0.5f , 0.8f } );
    flexnav.setupOutputMF( phit, phitACCEL, fis::trimf, (const real_t []){ 0.6f , 1.0f , 1.4f } );
    flexnav.setupOutputMF( thetat, thetatGYRO, fis::trimf, (const real_t []){ -0.4f ,0.0f , 0.4f } );
    flexnav.setupOutputMF( thetat, thetatBOTH, fis::trimf, (const real_t []){ 0.2f , 0.5f, 0.8f } );
    flexnav.setupOutputMF( thetat, thetatACCEL, fis::trimf, (const real_t []){ 0.6f , 1.0f , 1.4f } );

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
    cout << ffmath::mod( 3.1f, 2.5f ) << endl;
    cout << ffmath::rem( 3.1f, 2.5f ) << endl;
    cout << ffmath::wrapToPi( -4.5 ) << endl;
    cout << ffmath::getNan() << endl;
    cout << ffmath::getInf() << endl;
    cout << ffmath::tgamma( 2.5f ) << endl;
    cout << ffmath::lgamma( 2.5f ) << endl;

    cout << ffmath::assoc_laguerre( 1, 10, 0.5f ) << endl;
    cout << ffmath::assoc_laguerre( 2, 10, 0.5f ) << endl;


    cout << ffmath::assoc_legendre( 2, 0, 0.5f ) << endl;
    cout << ffmath::assoc_legendre( 2, 1, 0.5f ) << endl;
    cout << ffmath::assoc_legendre( 2, 2, 0.5f ) << endl;


    std::cout << ffmath::beta(0.1f, 0.2f) << std::endl;
    std::cout << "comp_ellint_1" << std::endl;
    std::cout << ffmath::comp_ellint_1( 0.0f ) << std::endl;
    std::cout << ffmath::comp_ellint_1( 0.5f ) << std::endl;
    std::cout << "comp_ellint_2" << std::endl;
    std::cout << ffmath::comp_ellint_2( 0.0f ) << std::endl;
    std::cout << ffmath::comp_ellint_2( 1.0f ) << std::endl;
    std::cout << ffmath::comp_ellint_2( 0.5f ) << std::endl;
    std::cout << "comp_ellint_3" << std::endl;
    std::cout << ffmath::comp_ellint_3( 0.5f, 0.0f ) << std::endl;
    std::cout << ffmath::comp_ellint_3( 0.0f, 0.0f ) << std::endl;
    std::cout << ffmath::comp_ellint_3( 0.5, 1.0f ) << std::endl;
    std::cout << "ellint_2" << std::endl;
    std::cout << ffmath::ellint_1( 0.0f, ffmath::FFP_PI_2 ) << std::endl;
    std::cout << ffmath::ellint_1( 0.0f, -ffmath::FFP_PI_2 ) << std::endl;
    std::cout << ffmath::ellint_1( 0.7f, 0.0f ) << std::endl;
    std::cout << "ellint_2" << std::endl;
    std::cout << ffmath::ellint_2( 0.0f , ffmath::FFP_PI_2 ) << std::endl;
    std::cout << ffmath::ellint_2( 0.0f , -ffmath::FFP_PI_2 ) << std::endl;
    std::cout << ffmath::ellint_2( 0.7f, 0.0f ) << std::endl;
    std::cout << ffmath::ellint_2( 1.0f, ffmath::FFP_PI_2 ) << std::endl;
    std::cout << "ellint_3" << std::endl;
    std::cout << ffmath::ellint_3( 0.0f, 0.0f, ffmath::FFP_PI_2 ) << std::endl;

    std::cout << "midpoint" << std::endl;
    std::cout << ffmath::midpoint( 3.0f, 6.0f ) << std::endl;
    std::cout << ffmath::midpoint( 6.0f, 3.0f ) << std::endl;
    std::cout << ffmath::midpoint( 6.56f, 7.23f ) << std::endl;
    std::cout << ffmath::midpoint( 2.0f, 3.0f ) << std::endl;
    std::cout << "lerp" << std::endl;
    std::cout << ffmath::lerp( 5.0f, 10.0f, -2.0f ) << std::endl;
    std::cout << ffmath::lerp( 5.0f, 10.0f, -1.5f ) << std::endl;
    std::cout << ffmath::lerp( 5.0f, 10.0f, -1.0f ) << std::endl;
    std::cout << ffmath::lerp( 5.0f, 10.0f, -0.5f ) << std::endl;

    std::cout << "expint" << std::endl;
    std::cout << ffmath::expint( 0.0f ) << std::endl;
    std::cout << ffmath::expint( 1.0f ) << std::endl;
    std::cout << -ffmath::exp(1.0f)*ffmath::expint( -1.0f ) << std::endl;
    //for (float x{1.f}; x < 8.8f; x += 0.3565f){
    //    std::cout << ffmath::expint( x ) << std::endl;
    //}
    std::cout << "hermite" << std::endl;
    std::cout << ffmath::hermite( 3, 10 ) << std::endl;
    std::cout << ffmath::hermite( 4, 10 ) << std::endl;

    std::cout << "laguerre" << std::endl;
    std::cout << ffmath::laguerre( 1, 0.5 ) << std::endl;
    std::cout << ffmath::laguerre( 2, 0.5 ) << std::endl;
    std::cout << ffmath::laguerre( 3, 0 ) << std::endl;

    std::cout << "legendre" << std::endl;
    std::cout << ffmath::legendre( 3, 0.25 ) << std::endl;
    std::cout << ffmath::legendre( 4, 0.25 ) << std::endl;

    std::cout << "riemann_zeta" << std::endl;
    for (const float x : {-39.0f, -2.0f, -1.0f, 0.0f, 1.0f, 0.5f, 2.0f})
        std::cout << x << " -> "<< ffmath::riemann_zeta(x) << std::endl;

    std::cout << "sph_bessel" << std::endl;
    std::cout << ffmath::sph_bessel( 1, 1.2345f ) << std::endl;
    std::cout << "sph_neumann" << std::endl;
    std::cout << ffmath::sph_neumann( 1, 1.2345f ) << std::endl;

    std::cout << "cyl_bessel_i" << std::endl;
    std::cout << ffmath::cyl_bessel_i( 0.0f, 1.2345f ) << std::endl;
    std::cout << ffmath::cyl_bessel_i( 1.0f, 1.2345f ) << std::endl;

    std::cout << "cyl_bessel_k" << std::endl;
    std::cout << ffmath::cyl_bessel_k( 0.5f, 1.2345f ) << std::endl;

    std::cout << "cyl_bessel_j" << std::endl;
    std::cout << ffmath::cyl_bessel_j( 0.0f, 1.2345f ) << std::endl;
    std::cout << ffmath::cyl_bessel_j( 0.0f, 1100.0f  ) << std::endl;

    std::cout << "sph_legendre" << std::endl;
    std::cout << ffmath::sph_legendre( 3, 0, 2.2345f ) << std::endl;
    std::cout << "copysgn" << std::endl;
    std::cout << ffmath::copysign( 1.0f, -2.0f ) << std::endl;
    std::cout << ffmath::copysign( ffmath::getNan(), -2.0f ) << std::endl;
    std::cout << ffmath::copysign( ffmath::getInf(), -2.0f ) << std::endl;
}

void test_mat( void )
{
    cout<<"MAT TEST"<<endl;
    mat<2,2> x(
                1.0f, 2.0f,
                3.0f, 4.0f
              );
    mat<2,1> y(
                1.0f,
                3.0f
              );
    mat<2,1> z(
                2.0f,
                2.0f
              );

    mat<2,1> result = 4.5f - z - 2.0f*x*y - 3.0f;
    mat<4,4> I( MAT_IDENTITY );
    x*=x;
    result.display();
    x.display();
    cout<< x(2) << endl;
    auto j = !y;
    j.display();
    I.display();
    auto ix = x.inv();
    x.display();
    ix.display();
    (x*ix).display();

}

void test_interp1( void )
{
    real_t tx[] = { 1.0f, 2.0f, 3.0f, 4.0f };
    real_t ty[] = { 5.0f, 9.0f, 12.0f, 15.0f };
    interp1 interpolation( tx, ty );

    cout << "linear" << endl;
    interpolation.setMethod( INTERP1_LINEAR );
    cout << interpolation.get( 2.5f ) << endl;
    cout << interpolation.get( 3.1f ) << endl;
    cout << interpolation.get( 0.5f ) << endl;
    cout << interpolation.get( 5.0f ) << endl;

    cout << "sine" << endl;
    interpolation.setMethod( INTERP1_SINE );
    cout << interpolation.get( 2.5f ) << endl;
    cout << interpolation.get( 3.1f ) << endl;
    cout << interpolation.get( 0.5f ) << endl;
    cout << interpolation.get( 5.0f ) << endl;

    cout << "cubic" << endl;
    interpolation.setMethod( INTERP1_CUBIC );
    cout << interpolation.get( 2.5f ) << endl;
    cout << interpolation.get( 3.1f ) << endl;
    cout << interpolation.get( 0.5f ) << endl;
    cout << interpolation.get( 5.0f ) << endl;

    cout << "hermite" << endl;
    interpolation.setMethod( INTERP1_HERMITE );
    cout << interpolation.get( 2.5f ) << endl;
    cout << interpolation.get( 3.1f ) << endl;
    cout << interpolation.get( 0.5f ) << endl;
    cout << interpolation.get( 5.0f ) << endl;

    cout << "nearest" << endl;
    interpolation.setMethod( INTERP1_NEAREST);
    cout << interpolation.get( 2.5f ) << endl;
    cout << interpolation.get( 3.1f ) << endl;
    cout << interpolation.get( 0.5f ) << endl;
    cout << interpolation.get( 5.0f ) << endl;

    cout << "next" << endl;
    interpolation.setMethod( INTERP1_NEXT);
    cout << interpolation.get( 2.5f ) << endl;
    cout << interpolation.get( 3.1f ) << endl;
    cout << interpolation.get( 0.5f ) << endl;
    cout << interpolation.get( 5.0f ) << endl;

    cout << "previous" << endl;
    interpolation.setMethod( INTERP1_PREVIOUS);
    cout << interpolation.get( 2.5f ) << endl;
    cout << interpolation.get( 3.1f ) << endl;
    cout << interpolation.get( 0.5f ) << endl;
    cout << interpolation.get( 5.0f ) << endl;

    cout << "spline" << endl;
    interpolation.setMethod( INTERP1_SPLINE);
    cout << interpolation.get( 2.5f ) << endl;
    cout << interpolation.get( 3.1f ) << endl;
    cout << interpolation.get( 0.5f ) << endl;
    cout << interpolation.get( 5.0f ) << endl;

    cout << "cspline" << endl;
    interpolation.setMethod( INTERP1_CONSTRAINED_SPLINE);
    cout << interpolation.get( 2.5f ) << endl;
    cout << interpolation.get( 3.1f ) << endl;
    cout << interpolation.get( 0.5f ) << endl;
    cout << interpolation.get( 5.0f ) << endl;


    real_t xdat[] = { 1.0f, 6.0f, 11.0f, 16.0f, 21.0f, 26.0f, 31.0f, 36.0f };
    real_t ydat[] = { 59.6870f,  44.5622f, -0.8642f , 0.8725f, -2.3016f, -50.3095f, -54.5966f, 37.9036f };
    interp1 interpolation2( xdat, ydat );
    interpolation2.setMethod( INTERP1_HERMITE );
    cout << "spline test 2" << endl;
    for( int i = 1; i<36;i++) {

        cout << i << "  " <<  interpolation2.get( i ) << endl;
    }
}

struct thing{
    int a;
    float b;
};


bool operator<(const thing& lhs, const thing& rhs) {
    return ( lhs.a < rhs. a);
}

int main()
{


    thing things[] = {
        {1,0},
        {-2,5},
        {8,4},
        {-2, 2},
    };

    for ( auto i : things ) {
        std::cout <<  i.a << " , "<< i.b << std::endl;
    }
    std::cout << "sorting "<< std::endl;
    algorithm::sort( things );



    for ( auto i : things ) {
        std::cout <<  i.a << " , "<< i.b << std::endl;
    }

    std::cout << "reverse "<< std::endl;
    algorithm::reverse( things );
    for ( auto i : things ) {
        std::cout <<  i.a << " , "<< i.b << std::endl;
    }
    std::cout << "rotate "<< std::endl;
    algorithm::rotate( things );
    for ( auto i : things ) {
        std::cout <<  i.a << " , "<< i.b << std::endl;
    }

    return 0;
    test_crc();
    test_fp16();
    test_ltisys();
    test_tdl();

    test_mat();
    test_interp1();
    test_fis();
    test_fis2();
    test_fis3();
    test_ffmath();

    return 0;
}