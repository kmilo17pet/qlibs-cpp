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
        (void)setGains( kc, kc/ti, kc*td );
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
        nextGains = { kc, ki, kd };
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
        nextGains = g;
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
        nextGains = { Kc*tmp, Kc/( ti*tmp ), Kc*( td/tmp ) };
        Kc = nextGains.Kc;
        Ki = nextGains.Ki;
        Kd = nextGains.Kd;
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
bool pidAutoTuning::getTimeConstant( real_t& tau,
                                     const real_t abs_z,
                                     const real_t dt )
{
    const bool stable = ( abs_z < 1.0_re );

    if ( stable && ( abs_z > 1e-12 ) ) {
        tau = -dt/ffmath::log( abs_z );
    }

    return stable;
}
/*============================================================================*/
void pidAutoTuning::computeParamsN1( systemParams& p,
                                     const estimationParams& e,
                                     const real_t dt )
{
    const real_t b1 = e.theta[ 0 ];
    const real_t a1 = e.theta[ 1 ];
    const real_t abs_z1 = ffmath::absf( a1 );

    p.gain = b1/( 1.0_re + a1 );
    p.tau2 = 0.0_re;
    p.stable = getTimeConstant( p.tau1, abs_z1, dt );
}
/*============================================================================*/
void pidAutoTuning::computeParamsN2( systemParams& p,
                                     const estimationParams& e,
                                     const real_t dt )
{
    const real_t b1 = e.theta[ 0 ];
    const real_t b2 = e.theta[ 1 ];
    const real_t a1 = e.theta[ 2 ];
    const real_t a2 = e.theta[ 3 ];

    const real_t a1a1 = a1*a1;
    const real_t a2_4 = 4.0_re*a2;
    const real_t d = a1a1 - a2_4;

    p.gain = ( b1 + b2 )/( 1.0_re + a1 + a2 );
    p.stable = true;

    if ( d >= 0.0_re ) {
        const real_t sqrt_d = ffmath::sqrt( d );
        const real_t  z1 = (-a1 + sqrt_d)*0.5_re;
        const real_t  z2 = (-a1 - sqrt_d)*0.5_re;
        const real_t  abs_z1 = ffmath::absf( z1 );
        const real_t  abs_z2 = ffmath::absf( z2 );

        p.stable = getTimeConstant( p.tau1, abs_z1, dt );
        p.stable &= getTimeConstant( p.tau2, abs_z2, dt );
    }
    else {
        const real_t rPart = -0.5_re*a1;
        const real_t iPart = -0.5_re*ffmath::sqrt( a2_4 - a1a1 );
        const real_t r = ffmath::hypot( rPart, iPart );
        p.tau2 = 0.0_re;
        p.stable = getTimeConstant( p.tau1, r, dt );
    }
}
/*============================================================================*/
bool pidAutoTuning::resetEstimation( estimationParams &p,
                                     const real_t r )
{
    bool rst = ( r < 0.0_re ) || ( r > 10000.0_re );
    if ( rst ) {
        constexpr real_t r_alfa = 0.1_re;
        p.P[0][0] = r_alfa; p.P[0][1] = 0.0_re; p.P[0][2] = 0.0_re; p.P[0][3] = 0.0_re;
        p.P[1][0] = 0.0_re; p.P[1][1] = r_alfa; p.P[1][2] = 0.0_re; p.P[1][3] = 0.0_re;
        p.P[2][0] = 0.0_re; p.P[2][1] = 0.0_re; p.P[2][2] = r_alfa; p.P[2][3] = 0.0_re;
        p.P[3][0] = 0.0_re; p.P[3][1] = 0.0_re; p.P[3][2] = 0.0_re; p.P[3][3] = r_alfa;
    }
    return rst;
}
/*============================================================================*/
void pidAutoTuning::estimationStepN1( estimationParams &p,
                                      const real_t y )
{
    constexpr real_t reg = 0.0001_re;
    const real_t f = 1.0_re/p.lambda;
    const real_t phi[2] = { p.u1, p.y2 };
    const real_t P_phi[2] = {
        ( p.P[0][0]*phi[0] ) + ( p.P[0][1]*phi[1] ),
        ( p.P[1][0]*phi[0] ) + ( p.P[1][1]*phi[1] ),
    };
    const real_t r = p.lambda + ( ( phi[0]*P_phi[0] ) + ( phi[1]*P_phi[1] ) );
    if ( resetEstimation( p, r ) ) {
        return;
    }
    const real_t inv_r = 1.0_re/r;
    const real_t L[4] ={ P_phi[0]*inv_r, P_phi[1]*inv_r };
    const real_t e  = y - ( ( phi[0]*p.theta[0] ) + ( phi[1]*p.theta[1] ) );
    p.theta[0] += L[0]*e;
    p.theta[1] += L[1]*e;
    p.P[0][0] = f*( p.P[0][0] - ( L[0]*P_phi[0] ) ) + reg;
    p.P[0][1] = f*( p.P[0][1] - ( L[0]*P_phi[1] ) );
    p.P[1][0] = f*( p.P[1][0] - ( L[1]*P_phi[0] ) );
    p.P[1][1] = f*( p.P[1][1] - ( L[1]*P_phi[1] ) ) + reg;
}
/*============================================================================*/
void pidAutoTuning::estimationStepN2( estimationParams &p,
                                      const real_t y )
{
    constexpr real_t reg = 0.0001_re;
    const real_t f = 1.0_re/p.lambda;
    const real_t phi[4] = { p.u1, p.u2, p.y1, p.y2 };
    const real_t P_phi[4] = {
         ( p.P[0][0]*phi[0] ) + ( p.P[0][1]*phi[1] ) + ( p.P[0][2]*phi[2] ) + ( p.P[0][3]*phi[3] ),
         ( p.P[0][1]*phi[0] ) + ( p.P[1][1]*phi[1] ) + ( p.P[1][2]*phi[2] ) + ( p.P[1][3]*phi[3] ),
         ( p.P[0][2]*phi[0] ) + ( p.P[1][2]*phi[1] ) + ( p.P[2][2]*phi[2] ) + ( p.P[2][3]*phi[3] ),
         ( p.P[0][3]*phi[0] ) + ( p.P[1][3]*phi[1] ) + ( p.P[2][3]*phi[2] ) + ( p.P[3][3]*phi[3] ),
    };
    const real_t r = p.lambda + ( ( phi[0]*P_phi[0] ) + ( phi[1]*P_phi[1] ) + ( phi[2]*P_phi[2] ) + ( phi[3]*P_phi[3]) );
    if ( resetEstimation( p, r ) ) {
        return;
    }
    const real_t inv_r = 1.0_re/r;
    const real_t L[4] ={ P_phi[0]*inv_r, P_phi[1]*inv_r, P_phi[2]*inv_r, P_phi[3]*inv_r };
    const real_t e  = y - ( ( phi[0]*p.theta[0] ) + ( phi[1]*p.theta[1] ) + ( phi[2]*p.theta[2] ) + ( phi[3]*p.theta[3] ) );
    p.theta[0] += L[0]*e;
    p.theta[1] += L[1]*e;
    p.theta[2] += L[2]*e;
    p.theta[3] += L[3]*e;
    p.P[0][0] = f*( p.P[0][0] - ( L[0]*P_phi[0] ) ) + reg;
    p.P[0][1] = f*( p.P[0][1] - ( L[0]*P_phi[1] ) );
    p.P[0][2] = f*( p.P[0][2] - ( L[0]*P_phi[2] ) );
    p.P[0][3] = f*( p.P[0][3] - ( L[0]*P_phi[3] ) );
    p.P[1][1] = f*( p.P[1][1] - ( L[1]*P_phi[1] ) ) + reg;
    p.P[1][2] = f*( p.P[1][2] - ( L[1]*P_phi[2] ) );
    p.P[1][3] = f*( p.P[1][3] - ( L[1]*P_phi[3] ) );
    p.P[2][2] = f*( p.P[2][2] - ( L[2]*P_phi[2] ) ) + reg;
    p.P[2][3] = f*( p.P[2][3] - ( L[2]*P_phi[3] ) );
    p.P[3][3] = f*( p.P[3][3] - ( L[3]*P_phi[3] ) ) + reg;
}
/*============================================================================*/
bool pidAutoTuning::step( const real_t u,
                          const real_t y,
                          const real_t dt ) noexcept
{
    bool ready = false;

    estimationStep( estParams, y );
    estParams.y2 = estParams.y1;
    estParams.y1 = -y;
    estParams.u2 = estParams.u1;
    estParams.u1 = u;
    computeParameters( sysParams, estParams, dt );
    if ( sysParams.stable && ( it > 0UL ) ) {
        if ( ( 0UL == --it ) && ( pidAutoTuning::UNDEFINED != it ) ) {
            ready = true;
        }
    }
    return ready;
}
/*============================================================================*/
bool pidAutoTuning::enable2ndOrderEstimates( const bool en ) noexcept
{
    if ( en ) {
        estimationStep = estimationStepN2;
        computeParameters = computeParamsN2;
    }
    else {
        estimationStep = estimationStepN1;
        computeParameters = computeParamsN1;
    }
}
/*============================================================================*/
pidGains pidAutoTuning::getEstimates( void ) const noexcept
{
    pidGains gains = { 0.0_re, 0.0_re, 0.0_re };
    const real_t tau = sysParams.tau1 + sysParams.tau2;
    const real_t k = sysParams.gain;
    const real_t td = tau*0.1_re;

    switch ( type ) {
        case pidType::PID_TYPE_P:
            gains.Kc = speed*( 1.03_re/k )*( ( tau/td ) + 0.34_re );
            break;
        case pidType::PID_TYPE_PD:
            gains.Kc = speed*( 1.24_re/k )*( ( tau/td ) + 0.129_re );
            gains.Kd = gains.Kc*( 0.27_re*td )*( tau - 0.324_re*td )/( tau + 0.129_re*td );
            break;
        case pidType::PID_TYPE_PI:
            gains.Kc = speed*( 0.9_re/k )*( ( tau/td ) + 0.092_re );
            gains.Ki = gains.Kc*( tau + 2.22_re*td )/( 3.33_re*td*( tau + 0.092_re*td ) );
            break;
        case pidType::PID_TYPE_PID:
            gains.Kc = speed*( 1.35_re/k )*( ( tau/td ) + 0.185_re );
            gains.Ki = gains.Kc*( tau + 0.611_re*td )/( 2.5_re*td*( tau + 0.185_re*td ) );
            gains.Kd = ( 0.37_re*gains.Kc*td*tau )/( tau + 0.185_re*td );
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
    real_t th0 = 0.25_re;   // b0
    real_t th1 = -0.077_re; // b1
    real_t th2 = -1.1_re;   // a0
    real_t th3 = 0.28_re;   // b1
    if ( &pidAutoTuning::computeParamsN1 == computeParameters ) {
        real_t k_i, T_i;
        k_i = current.Kc/0.9_re;
        T_i = ( 0.27_re*k_i )/current.Ki;
        /*cstat -CERT-FLP32-C_b*/
        th1 = -ffmath::exp( -dt/T_i ); //a1
        /*cstat +CERT-FLP32-C_b*/
        th0 = k_i*( 1.0_re + th1 );  //b1
    }

    estParams = {
        { th0, th1, th2, th3 },
        {
            { 1000.0_re, 0.0_re,    0.0_re,    0.0_re    },
            { 0.0_re,    1000.0_re, 0.0_re,    0.0_re    },
            { 0.0_re,    0.0_re,    1000.0_re, 0.0_re    },
            { 0.0_re,    0.0_re,    0.0_re,    1000.0_re },
        },
        0.9898_re,
        0.0_re, 0.0_re, 0.0_re, 0.0_re
    };
    sysParams = {
        0.0_re,
        0.0_re,
        0.0_re,
        false,
    };
    it = 100UL;
    mu = 0.95_re;


    speed = 0.25_re;
    it = pidAutoTuning::UNDEFINED;
}
/*============================================================================*/
void pidController::adaptGains( const real_t u,
                                const real_t y ) noexcept
{
    if ( adapt->step( u, y, dt ) ) {
        nextGains = adapt->getEstimates();
    }
    Kc += gainBlend*( nextGains.Kc - Kc );
    Ki += gainBlend*( nextGains.Ki - Ki );
    Kd += gainBlend*( nextGains.Kd - Kd );
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
                                             const real_t lambda,
                                             const real_t Tb ) noexcept
{
    bool retValue = false;

    if ( nullptr != adapt ) {
        if ( adapt->isValidParam( Mu ) && adapt->isValidParam( Alpha ) && adapt->isValidParam( lambda ) && ( Tb > 0.0_re )  ) { //skipcq : CXX-C2022
            adapt->setMemoryFactor( lambda );
            adapt->setMomentum( Mu );
            adapt->setEstimatedControllerSpeed( Alpha );
            gainBlend = dt/Tb;
            if ( gainBlend > 1.0_re ) {
                gainBlend = 1.0_re;
            }
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
    (void)ffmath::inRangeCoerce( x, vMin, vMax );
    return x;
}
/*============================================================================*/
