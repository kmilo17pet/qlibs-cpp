#include <include/fis.hpp>

using namespace qlibs;

const size_t fis::instance::INFERENCE_ERROR = 0U;
/*============================================================================*/
bool fis::instance::setParameter( const fis::parameter p,
                                  const fis::paramValue x ) noexcept
{
    bool retValue = false;
    static const methods_fcn method[ 5 ] = {
        &Min, &Prod, &Max, &ProbOr, &Sum
    };
    switch ( p ) {
        case fis::parameter::FIS_Aggregation:
            if ( ( x >= fis::paramValue::FIS_MAX ) && ( x <= fis::paramValue::FIS_SUM ) ) {
                aggregate = method[ x ];
                retValue = true;
            }
            break;
        case fis::parameter::FIS_Implication:
            if ( x <= fis::paramValue::FIS_PROBOR ) {
                implicate = method[ x ];
                retValue = true;
            }
            break;
        case fis::parameter::FIS_AND:
            if ( x <= fis::paramValue::FIS_PROD ) {
                andOp = method[ x ];
                retValue = true;
            }
            break;
        case fis::parameter::FIS_OR:
            if ( ( x >= fis::paramValue::FIS_MAX ) && ( x <= fis::paramValue::FIS_PROBOR ) ) {
                orOp = method[ x ];
                retValue = true;
            }
            break;
        default:
            break;
    }

    return retValue;
}
/*============================================================================*/
bool fis::instance::setDeFuzzMethod( fis::deFuzzMethod m ) noexcept
{
    bool retValue = false;

    static const fis::deFuzzFunction method[ FIS_NUM_DEFUZZ ] = { &deFuzzCentroid,
                                                                  &deFuzzBisector,
                                                                  &deFuzzMOM,
                                                                  &deFuzzLOM,
                                                                  &deFuzzSOM,
                                                                  &deFuzzWtAverage,
                                                                  &deFuzzWtSum
                                                                };

    if ( m < FIS_NUM_DEFUZZ ) {
        /*cppcheck-suppress knownConditionTrueFalse */
        if ( ( ( Mamdani == xType ) && ( m <= som ) ) ||
             ( ( Sugeno == xType ) && ( m >= wtaver ) && ( m <= wtsum ) ) ||
             ( ( Tsukamoto == xType ) && ( wtaver == m ) ) ) {
            deFuzz = method[ m ];
            retValue = true;
        }
    }
    return retValue;
}
/*============================================================================*/
bool fis::instance::setup( const fis::type t,
                           fis::input * const inputs,
                           const size_t ni,
                           fis::output * const outputs,
                           const size_t no,
                           fis::mf * const mf_inputs,
                           const size_t nmi,
                           fis::mf * const mf_outputs,
                           const size_t nmo,
                           const fis::rules * const r,
                           const size_t n,
                           real_t *rWeights ) noexcept
{
    bool retValue = false;

    if ( ( t <= Tsukamoto ) &&
         ( nullptr != inputs ) && ( ni > 0U ) &&
         ( nullptr != outputs) && ( no > 0U ) &&
         ( nullptr != mf_inputs) && ( nmi > 0U ) &&
         ( nullptr != mf_outputs) && ( nmo > 0U ) &&
         ( nullptr != r ) && ( n > 0U )
       ) {
        xType = t;
        xInput = inputs;
        xOutput = outputs;
        nInputs = ni;
        nOutputs = no;
        inMF = mf_inputs;
        nMFInputs = nmi;
        outMF = mf_outputs;
        nMFOutputs = nmo;
        xRules = r;
        wi = rWeights;
        nRules = n;
        deFuzz = ( Mamdani == xType )? &deFuzzCentroid : &deFuzzWtAverage;
        (void)setParameter( FIS_AND, FIS_MIN );
        (void)setParameter( FIS_OR, FIS_MAX );
        (void)setParameter( FIS_Implication, FIS_MIN );
        (void)setParameter( FIS_Aggregation, FIS_MAX );
        for ( size_t i = 0U ; i < nOutputs ; ++i ) {
            xOutput[ i ].owner = this;
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool fis::instance::setupInput( const fis::tag t,
                                const real_t Min,
                                const real_t Max ) noexcept
{
    bool retVal = false;

    if ( ( nullptr != xInput ) && ( t >= 0 ) ) {
        /*cstat -CERT-STR34-C*/
        xInput[ t ].min = Min;
        xInput[ t ].max = Max;
        /*cstat +CERT-STR34-C*/
        retVal = true;
    }

    return retVal;
}
/*============================================================================*/
bool fis::instance::setupOutput( const fis::tag t,
                                 const real_t Min,
                                 const real_t Max ) noexcept
{
    bool retVal = false;

    if ( ( nullptr != xInput ) && ( t >= 0 ) ) {
        /*cstat -CERT-STR34-C*/
        xOutput[ t ].min = Min;
        xOutput[ t ].max = Max;
        /*cstat -CERT-FLP36-C*/
        xOutput[ t ].res = ( xOutput[ t ].max - xOutput[ t ].min )/static_cast<real_t>( nPoints );
        /*cstat +CERT-STR34-C +CERT-FLP36-C*/
        retVal = true;
    }

    return retVal;
}
/*============================================================================*/
bool fis::output::storeAggregatedRegion( real_t *xData,
                                         real_t *yData,
                                         const size_t n ) noexcept
{
    bool retVal = false;
    if ( nullptr != owner ) {
        if ( ( nullptr != xData ) && ( nullptr != yData ) && ( n >= owner->getNumberOfPoints() ) ) {
            xag = xData;
            yag = yData;
            retVal = true;
        }
    }
    return retVal;
}
/*============================================================================*/
bool fis::instance::setInput( const fis::tag t,
                              const real_t value ) noexcept
{
    bool retVal = false;

    if ( ( nullptr != xInput ) && ( t >= 0 ) ) {
        /*cstat -CERT-STR34-C*/
        xInput[ t ].value = value;
        /*cstat +CERT-STR34-C*/
        retVal = true;
    }

    return retVal;
}
/*============================================================================*/
bool fis::instance::getOutput( const fis::tag t,
                               real_t &value ) const noexcept
{
    bool retVal = false;

    if ( ( nullptr != xOutput ) && ( t >= 0 ) ) {
        /*cstat -CERT-STR34-C*/
        value = xOutput[ t ].value;
        /*cstat +CERT-STR34-C*/
        retVal = true;
    }

    return retVal;
}
/*============================================================================*/
bool fis::instance::setMF( fis::mf *m,
                           const fis::tag io,
                           const fis::tag mf,
                           const fis::shapeMF s,
                           fis::mfFunction customMf,
                           const real_t *cp,
                           const real_t h ) noexcept
{
    bool retValue = false;
    static const fis::mfFunction fShape[ FIS_NUM_MFS ] = { &ConstantMF,
        /* Conventional membership functions, applies on any antecedent*/
        &TriMF, &TrapMF, &GBellMF, &GaussMF, &Gauss2MF, &SigMF, &DSigMF, &PSigMF,
        &PiMF, &SMF, &ZMF, &SingletonMF, &ConcaveMF, &SpikeMF, &LinSMF, &LinZMF,
        &RectangleMF, &CosineMF,
        /* Only for Sugeno consequents*/
        &ConstantMF, &LinearMF,
        /* Only for Tsukamoto consequents*/
        &TLinSMF, &TLinZMF, &TConcaveMF, &TSigMF, &TSMF, &TZMF, &TRampMF
    };

    if ( ( io >= 0 ) && ( mf >= 0 ) && ( s < FIS_NUM_MFS ) ) {
        /*cstat -CERT-STR34-C*/
        m[ mf ].shape = ( nullptr != customMf ) ? customMf : fShape[ s ];
        m[ mf ].index = static_cast<size_t>( io );
        m[ mf ].points = cp;
        m[ mf ].fx = 0.0_re;
        m[ mf ].h = bound( h );
        /*cstat +CERT-STR34-C*/
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool fis::instance::fuzzify( void ) noexcept
{
    bool retValue = false;

    if ( ( nullptr != xInput ) && ( nullptr != inMF ) ) {
        /* truncate inputs to its own range */
        for ( size_t i = 0U ; i < nInputs ; ++i ) {
            xInput[ i ].value = bound( xInput[ i ].value, xInput[ i ].min, xInput[ i ].max );
        }
        /* evaluate input membership functions */
        for ( size_t i = 0U ; i < nMFInputs ; ++i ) {
            const size_t mfIndex = inMF[ i ].getIndex();
            (void)inMF[ i ].membership( &xInput[ mfIndex ] );
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
real_t fis::instance::parseFuzzValue( fis::mf * const mfIO,
                                      fis::rules index ) noexcept
{
    const bool neg = ( index < 0 );
    real_t y;

    if ( neg ) {
        index = -index;
    }
    /*cstat -CERT-INT30-C_a -CERT-STR34-C -CERT-INT32-C_a*/
    y = bound( mfIO[ index - 1 ].fx );
    /*cstat +CERT-INT30-C_a +CERT-STR34-C +CERT-INT32-C_a*/
    y = ( neg ) ? ( 1.0_re - y ) : y;

    return y;
}
/*============================================================================*/
fis::fuzzyOperator fis::instance::getFuzzOperator( void ) noexcept
{
    fis::fuzzyOperator oper;

    switch ( lastConnector ) {
        case Q_FIS_AND:
            oper = andOp;
            break;
        case Q_FIS_OR:
            oper = orOp;
            break;
        default:
            oper = &Sum;
            break;
    }

    return oper;
}
/*============================================================================*/
size_t fis::instance::inferenceAntecedent( size_t i ) noexcept
{
    fis::rules inIndex, MFInIndex, connector;
    fis::fuzzyOperator op;

    inIndex = xRules[ i ];
    /*cstat -CERT-INT30-C_a*/
    MFInIndex = xRules[ i + 1U ];
    connector = xRules[ i + 2U ];
    /*cstat +CERT-INT30-C_a*/
    op = getFuzzOperator();
    rStrength = op( rStrength, parseFuzzValue( inMF, MFInIndex ) );

    if ( ( inIndex < 0 ) || ( static_cast<size_t>( inIndex ) > nInputs ) ) {
        i = INFERENCE_ERROR;
    }
    else {
        if ( ( Q_FIS_AND == connector ) || ( Q_FIS_OR == connector ) ) {
            lastConnector = connector;
            inferenceState = &fis::instance::inferenceAntecedent;
            i += 2U;
        }
        else if ( Q_FIS_THEN == connector ) {
            inferenceState = &fis::instance::inferenceReachEnd;
            i += 2U;
        }
        else {
            i = INFERENCE_ERROR;
        }
    }

    return i;
}
/*============================================================================*/
size_t fis::instance::inferenceReachEnd( size_t i ) noexcept
{
    fis::rules connector;
    /*cstat -CERT-INT30-C_a*/
    connector = ( nOutputs > 1U )? xRules[ i + 2U ] : -1;
    /*cstat +CERT-INT30-C_a*/
    i += 2U;
    if ( Q_FIS_AND != connector ) {
        inferenceState = &fis::instance::inferenceAntecedent;
        lastConnector = -1;
        wi[ ruleCount ] = rStrength;
        if ( nullptr != ruleWeight ) {
            wi[ ruleCount ] *= bound( ruleWeight[ ruleCount ] );
        }
        rStrength = 0.0_re;
        ++ruleCount;
        --i;
    }

    return i;
}
/*============================================================================*/
size_t fis::instance::aggregationFindConsequent( size_t i ) noexcept
{
    while ( Q_FIS_THEN != xRules[ i++ ] ) {}
    aggregationState = &fis::instance::inferenceConsequent;
    --i;

    return i;
}
/*============================================================================*/
size_t fis::instance::inferenceConsequent( size_t i ) noexcept
{
    fis::rules outIndex, MFOutIndex, connector;
    bool neg = false;

    outIndex = xRules[ i ];
    /*cstat -CERT-INT30-C_a*/
    MFOutIndex = xRules[ i + 1U ];
    connector = ( nOutputs > 1U )? xRules[ i + 2U ] : -1;
    /*cstat +CERT-INT30-C_a*/
    if ( MFOutIndex < 0 ) {
        MFOutIndex = -MFOutIndex;
        neg = true;
    }
    MFOutIndex -= 1;

    if ( wi[ ruleCount ] > 0.0_re ) {
        /*cstat -CERT-STR34-C*/
        fis::output &o = xOutput[ outIndex ];
        fis::mf &m = outMF[ MFOutIndex ];
        /*cstat +CERT-STR34-C*/

        if ( Mamdani == xType ) {
            real_t v;
            v = m.membership( &o );
            v = ( neg )? ( 1.0_re - v ) : v;
            o.y = aggregate( o.y, implicate( wi[ ruleCount ], v ) );
        }
        else { /* Sugeno and Tsukamoto*/
            real_t zi;
            zi = m.membership( xInput, nInputs );
            o.v[ sum_wz ] += zi*wi[ ruleCount ];
            o.v[ sum_w ] += wi[ ruleCount ];
        }
    }

    i += 2U;
    if ( Q_FIS_AND != connector ) {
        aggregationState = &fis::instance::aggregationFindConsequent;
        ++ruleCount;
        --i;
    }

    return i;
}
/*============================================================================*/
void fis::instance::fuzzyAggregate( void ) noexcept
{
    if ( Q_FIS_RULES_BEGIN == xRules[ 0 ] ) {
        size_t i = 1U;

        aggregationState = &fis::instance::aggregationFindConsequent;
        ruleCount = 0U;
        while ( ( Q_FIS_RULES_END != xRules[ i ] ) && ( ruleCount < nRules ) ) {
            i = (this->*aggregationState)( i );
            if ( INFERENCE_ERROR == i ) {
                break;
            }
            ++i;
        }
    }
}
/*============================================================================*/
bool fis::instance::deFuzzify( void ) noexcept
{
    bool retValue = false;

    if ( ( nullptr != xOutput ) && ( nullptr != outMF ) ) {
        size_t i;

        for ( i = 0U ; i < nOutputs ; ++i ) {
            deFuzz( &xOutput[ i ] , FIS_DEFUZZ_INIT );
        }
        if ( Mamdani == xType ) {
            size_t k;

            for ( k = 0U ; k < nPoints ; ++k ) {
                for ( i = 0U ; i < nOutputs ; ++i ) { /* initialize*/
                    xOutput[ i ].y = 0.0_re;
                    xOutput[ i ].getNextX( k );
                    xOutput[ i ].value = xOutput[ i ].x;
                }
                fuzzyAggregate();
                for ( i = 0U ; i < nOutputs ; ++i ) {
                    deFuzz( &xOutput[ i ] , FIS_DEFUZZ_COMPUTE );
                    if ( nullptr != xOutput[ i ].xag ) { /*store aggregated*/
                        xOutput[ i ].xag[ k ] = xOutput[ i ].x;
                        xOutput[ i ].yag[ k ] = xOutput[ i ].y;
                    }
                }
            }
        }
        else { /*Sugeno and Tsukamoto systems*/
            for ( i = 0U ; i < nOutputs ; ++i ) { /* initialize*/
                xOutput[ i ].v[ sum_wz ] = 0.0_re; /*store sum wi*/
                xOutput[ i ].v[ sum_w ] = 0.0_re; /*store sum zi*wi*/
            }
            fuzzyAggregate();
            for ( i = 0U ; i < nOutputs ; ++i ) { /* initialize*/
                deFuzz( &xOutput[ i ] , FIS_DEFUZZ_COMPUTE );
            }
        }
        for ( i = 0U ; i < nOutputs ; ++i ) {
            xOutput[ i ].value = deFuzz( &xOutput[ i ] , FIS_DEFUZZ_END );
            xOutput[ i ].value = bound( xOutput[ i ].value, xOutput[ i ].min, xOutput[ i ].max );
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool fis::instance::setRuleWeights( real_t *rWeights ) noexcept
{
    bool retValue = false;
    if ( nullptr != rWeights ) {
        ruleWeight = rWeights;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool fis::instance::inference( void ) noexcept
{
    bool retValue = false;

    if ( nullptr != xRules  ) {
        size_t i = 0U;

        if ( Q_FIS_RULES_BEGIN == xRules[ 0 ] ) {
            inferenceState = &fis::instance::inferenceAntecedent;
            rStrength = 0.0_re;
            lastConnector = -1;
            ruleCount = 0U;
            i = 1U;
            while ( ( Q_FIS_RULES_END != xRules[ i ] ) && ( ruleCount < nRules ) ) {
                i = (this->*inferenceState)( i );
                if ( INFERENCE_ERROR == i ) {
                    break;
                }
                ++i;
            }
        }
        if ( ( Q_FIS_RULES_END == xRules[ i ] ) && ( ruleCount == nRules) ) {
            retValue = true;
        }
    }

    return retValue;
}
/*============================================================================*/
