#include <include/pid.hpp>
#include <include/ffmath.hpp>

using namespace qlibs;


const uint32_t pidAutoTuning::UNDEFINED = 0xFFFFFFFEUL;

/*============================================================================*/
bool pidController::setup( const real_t kc,
                           const real_t ki,
                           const real_t kd,
                           const real_t dT ) noexcept
{
    bool retValue = false;

    if ( dt > 0.0_re ) {
        dt = dT;
        isInitialized = true;
        (void)setDerivativeFilter( 0.98_re );
        (void)setEpsilon( REAL_MIN );
        (void)setGains( kc, ki, kd );
        (void)setSaturation( 0.0_re, 100.0_re );
        (void)setMode( pidMode::PID_AUTOMATIC );
        (void)setManualInput( 0.0_re );
        (void)setExtraGains( 1.0_re, 1.0_re );
        (void)setDirection( pidDirection::PID_FORWARD );
        (void)setReferenceWeighting( 1.0_re, 0.0_re );
        retValue = reset();
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setDirection( const pidDirection d ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        dir = d;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setParams( const real_t kc,
                               const real_t ti,
                               const real_t td ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        Kc = kc;
        Ki = kc/ti;
        Kd = kc*td;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setGains( const real_t kc,
                              const real_t ki,
                              const real_t kd ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        Kc = kc;
        Ki = ki;
        Kd = kd;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setGains( const pidGains &g ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        Kc = g.Kc;
        Ki = g.Ki;
        Kd = g.Kd;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setExtraGains( const real_t Kw,
                                   const real_t Kt ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        kw = Kw;
        kt = Kt;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setSaturation( const real_t Min,
                                   const real_t Max ) noexcept
{
    bool retValue = false;

    if ( isInitialized && ( Max > Min ) ) {
        sat_Min = Min;
        sat_Max = Max;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setSeries( void ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        real_t ti, td, tmp;

        ti = Kc/Ki;
        td = Kd/Kc;
        tmp = 1.0_re + ( td/ti );
        Kc = Kc*tmp;
        Ki = Kc/( ti*tmp );
        Kd = Kc*( td/tmp );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setEpsilon( const real_t eps ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        epsilon = eps;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setDerivativeFilter( const real_t Beta ) noexcept
{
    bool retValue = false;

    if ( isInitialized && ( Beta >= 0.0_re ) && ( Beta < 1.0_re ) ) {
        beta = Beta;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setDerivativeFilterTimeConstant( const real_t Tf ) noexcept
{
    bool retValue = false;

    if ( isInitialized && ( Tf >= 0.0_re ) ) {
        beta = ffmath::exp( -dt/Tf );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setMode( const pidMode Mode ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        mode = Mode;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setReferenceWeighting( const real_t gb,
                                           const real_t gc ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        b = saturate( gb, 0.0_re, 1.0_re );
        c = saturate( gc, 0.0_re, 1.0_re );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setManualInput( const real_t manualInput ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        mInput = manualInput;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::reset( void ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        init(); //internal controller state
        m_state.init();
        b_state.init();
        D = 0.0_re;
        u1 = 0.0_re;
        m = 0.0_re;
        uSat = 0.0_re;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setModelReferenceControl( const real_t &modelRef,
                                              const real_t Gamma,
                                              const real_t Alpha ) noexcept
{
    bool retValue = false;

    if ( isInitialized && ( Gamma > 0.0_re ) && ( Alpha > 0.0_re ) ) {
        m_state.init();
        alpha = Alpha;
        gamma = Gamma;
        yr = &modelRef;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::removeModelReferenceControl( void ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        yr = nullptr;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
real_t pidController::error( real_t w, real_t y, real_t k ) noexcept
{
    real_t e = ( k*w ) - y;

    if ( ffmath::absf( e ) <= epsilon ) {
        e = 0.0_re;
    }

    return e;
}
/*============================================================================*/
real_t pidController::control( const real_t w,
                               const real_t y ) noexcept
{
    real_t u = w;

    if ( isInitialized ) {
        real_t e, v, de, ie, bt, sw, kc, ki, kd;
        kc = Kc;
        ki = Ki;
        kd = Kd;
        if ( pidDirection::PID_BACKWARD == dir ) {
            kc = ( kc > 0.0_re ) ? -kc : kc;
            ki = ( ki > 0.0_re ) ? -ki : ki;
            kd = ( kd > 0.0_re ) ? -kd : kd;
        }
        e = error( w, y );
        de = derive( error( w, y, c ) , dt, false );
        ie = integrate( e + u1 , dt );
        D = de + beta*( D - de ); /*derivative filtering*/
        v  = ( kc*error( w, y, b ) ) + ( ki*ie ) + ( kd*D ); /*compute PID action*/
        if ( nullptr != yr ) {
            /*MRAC additive controller using the modified MIT rule*/
            real_t theta = 0.0_re;
            if ( ffmath::absf( u1 ) <= epsilon ) { /*additive anti-windup*/
                const real_t em = y - yr[ 0 ];
                const real_t delta = -gamma*em*yr[ 0 ]/( alpha + ( yr[ 0 ]*yr[ 0 ] ) );
                theta = m_state.integrate( delta /*+ c->u1*/, dt );
            }
            v += w*theta;
        }
        /*bumpless-transfer*/
        bt = ( kt*mInput ) + ( kw*( uSat - m ) );
        m = b_state.integrate( bt, dt );
        sw = ( pidMode::PID_AUTOMATIC == mode ) ? v : m;
        uSat = saturate( sw, sat_Min, sat_Max );
        u = uSat; /*output saturated*/
        u1 = kw*( u - v ); /*anti-windup feedback*/
        if ( nullptr != adapt ) {
            adaptGains( u, y );
        }
    }

    return u;
}
/*============================================================================*/
bool pidAutoTuning::step( const real_t u,
                          const real_t y,
                          const real_t dt ) noexcept
{
    real_t error , r, l0, l1;
    real_t lp00, lp01, lp10, lp11;
    real_t tmp1, tmp2;
    real_t gain, timeConstant;
    const real_t il = 1.0_re/l;
    bool ready = false;

    tmp1 = p00*uk;
    tmp2 = p11*yk;
    r = l + ( uk*( tmp1 - ( p10*yk ) ) ) - ( yk*( ( p01*uk ) - tmp2 ) );
    /*compute corrections*/
    l0 = ( tmp1 - ( p01*yk ) )/r;
    l1 = ( ( p10*uk ) - tmp2 )/r;
    error = y - ( ( b1*uk ) - ( a1*yk ) );
    /*fix estimations*/
    b1 += l0*error;
    a1 += l1*error;
    /*update covariances*/
    lp00 = il*p00;
    lp01 = il*p01;
    lp10 = il*p10;
    lp11 = il*p11;
    tmp1 = ( l0*uk ) - 1.0_re;
    tmp2 = ( l1*yk ) + 1.0_re;
    p00 = ( l0*lp10*yk ) - ( lp00*tmp1 ) + 1e-10_re;
    p01 = ( l0*lp11*yk ) - ( lp01*tmp1 );
    p10 = ( lp10*tmp2 ) - ( l1*lp00*uk );
    p11 = ( lp11*tmp2 ) - ( l1*lp01*uk ) + 1e-10_re;
    /*update I/O measurements*/
    yk = y;
    uk = u;
    gain = b1/( 1.0_re + a1 );
    timeConstant = -dt/ffmath::log( ffmath::absf( a1 ) );
    /*cstat -MISRAC++2008-5-14-1*/
    if ( isValidValue( timeConstant ) && isValidValue( gain ) && ( it > 0UL ) ) { /*no side effects here*/
    /*cstat +MISRAC++2008-5-14-1*/
        k = gain + ( mu*( k - gain ) );
        tao = timeConstant + ( mu*( tao - timeConstant ) );
        if ( ( 0UL == --it ) && ( pidAutoTuning::UNDEFINED != it ) ) {
            ready = true;
        }
    }

    return ready;
}
/*============================================================================*/
pidGains pidAutoTuning::getEstimates( void ) const noexcept
{
    pidGains gains = { 0.0_re, 0.0_re, 0.0_re };
    const real_t td = tao*0.1_re;

    switch ( type ) {
        case pidType::PID_TYPE_P:
            gains.Kc = speed*( 1.03_re/k )*( ( tao/td ) + 0.34_re );
            break;
        case pidType::PID_TYPE_PD:
            gains.Kc = speed*( 1.24_re/k )*( ( tao/td ) + 0.129_re );
            gains.Kd = gains.Kc*( 0.27_re*td )*( tao - 0.324_re*td )/( tao + 0.129_re*td );
            break;
        case pidType::PID_TYPE_PI:
            gains.Kc = speed*( 0.9_re/k )*( ( tao/td ) + 0.092_re );
            gains.Ki = gains.Kc*( tao + 2.22_re*td )/( 3.33_re*td*( tao + 0.092_re*td ) );
            break;
        case pidType::PID_TYPE_PID:
            gains.Kc = speed*( 1.35_re/k )*( ( tao/td ) + 0.185_re );
            gains.Ki = gains.Kc*( tao + 0.611_re*td )/( 2.5_re*td*( tao + 0.185_re*td ) );
            gains.Kd = ( 0.37_re*gains.Kc*td*tao )/( tao + 0.185_re*td );
            break;
        default:
            break;
    }

    return gains;
}
/*============================================================================*/
void pidAutoTuning::initialize( const pidGains current,
                                const real_t dt ) noexcept
{
    real_t k_i, T_i;

    l = 0.9898_re;
    p00 = 1000.0_re;
    p11 = 1000.0_re;
    p01 = 0.0_re;
    p10 = 0.0_re;
    uk  = 0.0_re;
    yk = 0.0_re;
    k = 0.0_re;
    tao = 0.0_re;
    it = 100UL;
    mu = 0.95_re;
    k_i = current.Kc/0.9_re;
    T_i = ( 0.27_re*k_i )/current.Ki;
    /*cstat -CERT-FLP32-C_b*/
    a1 = -ffmath::exp( -dt/T_i );
    /*cstat +CERT-FLP32-C_b*/
    b1 = k_i*( 1.0_re + a1 );
    speed = 0.25_re;
    it = pidAutoTuning::UNDEFINED;
}
/*============================================================================*/
void pidController::adaptGains( const real_t u,
                                const real_t y ) noexcept
{
    if ( adapt->step( u, y, dt ) ) {
        pidGains newGains = adapt->getEstimates();
        Kc = newGains.Kc;
        Ki = newGains.Ki;
        Kd = newGains.Kd;
    }
}
/*============================================================================*/
bool pidController::bindAutoTuning( pidAutoTuning &at ) noexcept
{
    bool retValue = false;

    if ( isInitialized ) {
        adapt = &at;
        adapt->initialize( { Kc, Ki, Kd }, dt );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::enableAutoTuning( const uint32_t tEnable ) noexcept
{
    bool retValue = false;

    if ( nullptr != adapt ) {
        adapt->enable( tEnable );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidController::isAutoTuningComplete( void ) const noexcept
{
    bool retValue = false;

    if ( nullptr != adapt ) {
        retValue = adapt->isComplete();
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setAutoTuningParameters( const real_t Mu,
                                             const real_t Alpha,
                                             const real_t lambda ) noexcept
{
    bool retValue = false;

    if ( nullptr != adapt ) {
        if ( adapt->isValidParam( Mu ) && adapt->isValidParam( Alpha ) && adapt->isValidParam( lambda ) ) { //skipcq : CXX-C2022
            adapt->setMemoryFactor( lambda );
            adapt->setMomentum( Mu );
            adapt->setEstimatedControllerSpeed( Alpha );
            retValue = true;
        }
    }

    return retValue;
}
/*============================================================================*/
bool pidController::setAutoTuningControllerType( const pidType t ) noexcept
{
    bool retValue = false;

    if ( nullptr != adapt ) {
        adapt->setEstimatedControllerType( t );
        //adapt->type = t;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool pidAutoTuning::isValidValue( const real_t x ) noexcept
{
     /*cstat -MISRAC++2008-5-14-1*/
     return ( !ffmath::isNan( x ) ) && ( x > 0.0_re ) && ( !ffmath::isInf( x ) );
     /*cstat +MISRAC++2008-5-14-1*/
}
/*============================================================================*/
real_t pidController::saturate( real_t x,
                                const real_t vMin,
                                const real_t vMax ) noexcept
{
    if ( x > vMax ) {
        x = vMax;
    }
    else if ( x < vMin ) {
        x = vMin;
    }
    else {
        /*nothing to do*/
    }

    return x;
}
/*============================================================================*/
