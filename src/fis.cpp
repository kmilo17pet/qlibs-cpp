#include "include/fis.hpp"
#include "include/mathex.hpp"


using namespace qlibs;

const size_t fis::INFERENCE_ERROR = 0U;
/*============================================================================*/
bool fis::setParameter( const fisParameter p, const fisParamValue x )
{
    bool retValue = false;
    static const methods_fcn method[ 5 ] = {
        &Min, &Prod, &Max, &ProbOr, &Sum
    };
    switch ( p ) {
        case fisParameter::FIS_Aggregation:
            if ( ( x >= fisParamValue::FIS_MAX ) && ( x <= fisParamValue::FIS_SUM ) ) {
                aggregate = method[ x ];
                retValue = true;
            }
            break;
        case fisParameter::FIS_Implication:
            if ( x <= fisParamValue::FIS_PROBOR ) {
                implicate = method[ x ];
                retValue = true;
            }
            break;
        case fisParameter::FIS_AND:
            if ( x <= fisParamValue::FIS_PROD ) {
                andOp = method[ x ];
                retValue = true;
            }
            break;
        case fisParameter::FIS_OR:
            if ( ( x >= fisParamValue::FIS_MAX ) && ( x <= fisParamValue::FIS_PROBOR ) ) {
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
bool fis::setDeFuzzMethod( fisDeFuzzMethod m )
{
    bool retValue = false;

    static const fisDeFuzzFunction method[ _NUM_DEFUZZ ] = { &deFuzzCentroid,
                                                            &deFuzzBisector,
                                                            &deFuzzMOM,
                                                            &deFuzzLOM,
                                                            &deFuzzSOM,
                                                            &deFuzzWtAverage,
                                                            &deFuzzWtSum
                                                          };

    if ( m < _NUM_DEFUZZ ) {
        if ( ( ( Mamdani == type ) && ( m <= som ) ) ||
             ( ( Sugeno == type ) && ( m >= wtaver ) && ( m <= wtsum ) ) ||
             ( ( Tsukamoto == type ) && ( wtaver == m ) )) {
            deFuzz = method[ m ];
            retValue = true;
        }
    }
    return retValue;
}
/*============================================================================*/
bool fis::setup( const fisType t,
                 fisInput * const inputs, const size_t ni,
                 fisOutput * const outputs, const size_t no,
                 fisMF * const mf_inputs, const size_t nmi,
                 fisMF * const mf_outputs, const size_t nmo,
                 const fisRules * const r, const size_t n,
                 real_t *rWeights )
{
    bool retValue = false;

    if ( ( t <= Tsukamoto ) && 
         ( nullptr != inputs ) && ( ni > 0U ) && 
         ( nullptr != outputs) && ( no > 0U ) && 
         ( nullptr != mf_inputs) && ( nmi > 0U ) &&
         ( nullptr != mf_outputs) && ( nmo > 0U ) &&
         ( nullptr != r ) && ( n > 0U )
       ) {
        type = t;
        input = inputs;
        output = outputs;
        nInputs = ni;
        nOutputs = no;
        inMF = mf_inputs;
        nMFInputs = nmi;
        outMF = mf_outputs;
        nMFOutputs = nmo;
        rules = r;
        wi = rWeights;
        nRules = n;
        deFuzz = ( Mamdani == type )? &deFuzzCentroid : &deFuzzWtAverage;
        (void)setParameter( FIS_AND, FIS_MIN );
        (void)setParameter( FIS_OR, FIS_MAX );
        (void)setParameter( FIS_Implication, FIS_MIN );
        (void)setParameter( FIS_Aggregation, FIS_MAX );
        for ( size_t i = 0U ; i < nOutputs ; ++i ) {
            //output[ i ].res = ( output[ i ].max - output[ i ].min )/static_cast<real_t>( nPoints );
            output[ i ].owner = this;
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool fis::setupInput( const fisTag t, const real_t Min, const real_t Max )
{
    bool retVal = false;

    if ( ( nullptr != input ) && ( t >= 0 ) ) {
        /*cstat -CERT-STR34-C*/
        input[ t ].min = Min;
        input[ t ].max = Max;
        /*cstat +CERT-STR34-C*/
        retVal = true;
    }

    return retVal;
}
/*============================================================================*/
bool fis::setupOutput( const fisTag t, const real_t Min, const real_t Max )
{
    bool retVal = false;

    if ( ( nullptr != input ) && ( t >= 0 ) ) {
        /*cstat -CERT-STR34-C*/
        output[ t ].min = Min;
        output[ t ].max = Max;
        output[ t ].res = ( output[ t ].max - output[ t ].min )/static_cast<real_t>( nPoints );
        /*cstat +CERT-STR34-C*/
        retVal = true;
    }

    return retVal;
}
/*============================================================================*/
bool fisOutput::storeAggregatedRegion( real_t *xData, real_t *yData, const size_t n )
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
bool fis::setInput( const fisTag t, const real_t value )
{
    bool retVal = false;

    if ( ( nullptr != input ) && ( t >= 0 ) ) {
        /*cstat -CERT-STR34-C*/
        input[ t ].value = value;
        /*cstat +CERT-STR34-C*/
        retVal = true;
    }

    return retVal;
}
/*============================================================================*/
bool fis::getOutput( const fisTag t, real_t &value ) const
{
    bool retVal = false;

    if ( ( nullptr != output ) && ( t >= 0 ) ) {
        /*cstat -CERT-STR34-C*/
        value = output[ t ].value;
        /*cstat +CERT-STR34-C*/
        retVal = true;
    }

    return retVal;
}
/*============================================================================*/
bool fis::setMF( fisMF *m, const fisTag io, const fisTag mf, const fisShapeMF s, fisMFFunction customMf, const real_t *cp, const real_t h )
{
    bool retValue = false;
    static const fisMFFunction fShape[ _NUM_MFS ] = { &ConstantMF,
        /* Conventional membership functions, applies on any antecedent*/
        &TriMF, &TrapMF, &GBellMF, &GaussMF, &Gauss2MF, &SigMF, &DSigMF, &PSigMF,
        &PiMF, &SMF, &ZMF, &SingletonMF, &ConcaveMF, &SpikeMF, &LinSMF, &LinZMF,
        &RectangleMF, &CosineMF,
        /* Only for Sugeno consequents*/
        &ConstantMF, &LinearMF,
        /* Only for Tsukamoto consequents*/
        &TLinSMF, &TLinZMF, &TConcaveMF, &TSigMF, &TSMF, &TZMF
    };

    if ( ( io >= 0 ) && ( mf >= 0 ) && ( s < _NUM_MFS ) ) {
        /*cstat -CERT-STR34-C*/
        m[ mf ].shape = ( nullptr != customMf ) ? customMf : fShape[ s ];
        m[ mf ].index = static_cast<size_t>( io );
        m[ mf ].points = cp;
        m[ mf ].fx = 0.0;
        m[ mf ].h = bound( h, 0.0, 1.0 );
        /*cstat +CERT-STR34-C*/
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
void fis::evalInputMFs( void )
{
    for ( size_t i = 0U ; i < nMFInputs ; ++i ) {
        inMF[ i ].evalMFAtIndex( input );
    }
}
/*============================================================================*/
void fis::truncateInputs( void )
{
    for ( size_t i = 0U ; i < nInputs ; ++i ) {
        input[ i ].value = bound( input[ i ].value, input[ i ].min, input[ i ].max );
    }
}
/*============================================================================*/
bool fis::fuzzify( void )
{
    bool retValue = false;

    if ( ( nullptr != input ) && ( nullptr != inMF ) ) {
        truncateInputs();
        evalInputMFs();
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
real_t fis::parseFuzzValue( fisMF * const mfIO, fisRules index )
{
    const bool neg = ( index < 0 );
    real_t y;

    if ( neg ) {
        index = -index;
    }
    /*cstat -CERT-INT30-C_a -CERT-STR34-C -CERT-INT32-C_a*/
    y = bound( mfIO[ index - 1 ].fx, 0.0, 1.0 );
    /*cstat +CERT-INT30-C_a +CERT-STR34-C +CERT-INT32-C_a*/
    y = ( neg ) ? ( 1.0 - y ) : y;

    return y;
}
/*============================================================================*/
fuzzyOperator fis::getFuzzOperator( void )
{
    fuzzyOperator oper;

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
size_t fis::inferenceAntecedent( size_t i )
{
    fisRules inIndex, MFInIndex, connector;
    fuzzyOperator op;

    inIndex = rules[ i ];
    /*cstat -CERT-INT30-C_a*/
    MFInIndex = rules[ i + 1U ];
    connector = rules[ i + 2U ];
    /*cstat +CERT-INT30-C_a*/
    op = getFuzzOperator();
    rStrength = op( rStrength, parseFuzzValue( inMF, MFInIndex ) );
 
    if ( ( inIndex < 0 ) || ( static_cast<size_t>( inIndex ) > nInputs ) ) {
        i = INFERENCE_ERROR;
    }
    else {
        if ( ( Q_FIS_AND == connector ) || ( Q_FIS_OR == connector ) ) {
            lastConnector = connector;
            inferenceState = &fis::inferenceAntecedent;
            i += 2U;
        }
        else if ( Q_FIS_THEN == connector ) {
            inferenceState = &fis::inferenceReachEnd;
            i += 2U;
        }
        else {
            i = INFERENCE_ERROR;
        }
    }

    return i;
}
/*============================================================================*/
size_t fis::inferenceReachEnd( size_t i )
{
    fisRules connector;
    /*cstat -CERT-INT30-C_a*/
    connector = ( nOutputs > 1U )? rules[ i + 2U ] : -1;
    /*cstat +CERT-INT30-C_a*/
    i += 2U;
    if ( Q_FIS_AND != connector ) {
        inferenceState = &fis::inferenceAntecedent;
        lastConnector = -1;
        wi[ ruleCount ] = rStrength;
        if ( nullptr != ruleWeight ) {
            wi[ ruleCount ] *= bound( ruleWeight[ ruleCount ], 0.0, 1.0 );
        }
        rStrength = 0.0;
        ++ruleCount;
        --i;
    }

    return i;
}
/*============================================================================*/
size_t fis::aggregationFindConsequent( size_t i )
{
    while ( Q_FIS_THEN != rules[ i++ ] ) {}
    aggregationState = &fis::inferenceConsequent;

    return --i;
}
/*============================================================================*/
size_t fis::inferenceConsequent( size_t i )
{
    fisRules outIndex, MFOutIndex, connector;
    bool neg = false;

    outIndex = rules[ i ];
    /*cstat -CERT-INT30-C_a*/
    MFOutIndex = rules[ i + 1U ];
    connector = ( nOutputs > 1U )? rules[ i + 2U ] : -1;
    /*cstat +CERT-INT30-C_a*/
    if ( MFOutIndex < 0 ) {
        MFOutIndex = -MFOutIndex;
        neg = true;
    }
    MFOutIndex -= 1;

    if ( wi[ ruleCount ] > 0.0 ) {
        /*cstat -CERT-STR34-C*/
        fisOutput &o = output[ outIndex ];
        fisMF &m = outMF[ MFOutIndex ];
        /*cstat +CERT-STR34-C*/
        
        if ( Mamdani == type ) {
            real_t v;
            v = m.evalMF( o );
            v = ( neg )? ( 1.0 - v ) : v;
            o.y = aggregate( o.y, implicate( wi[ ruleCount ], v ) );
        }
        else { /* Sugeno and Tsukamoto*/
            real_t zi;
            zi = m.evalMF( input, nInputs );
            o.v[ sum_wz ] += zi*wi[ ruleCount ];
            o.v[ sum_w ] += wi[ ruleCount ];
        }
    }

    i += 2U;
    if ( Q_FIS_AND != connector ) {
        aggregationState = &fis::aggregationFindConsequent;
        ++ruleCount;
        --i;
    }

    return i;
}
/*============================================================================*/
void fis::fuzzyAggregate( void )
{
    if ( Q_FIS_RULES_BEGIN == rules[ 0 ] ) {
        size_t i = 1U;

        aggregationState = &fis::aggregationFindConsequent;
        ruleCount = 0U;
        while ( ( Q_FIS_RULES_END != rules[ i ] ) && ( ruleCount < nRules ) ) {
            i = (this->*aggregationState)( i );
            if ( INFERENCE_ERROR == i ) {
                break;
            }
            ++i;
        }
    }
}
/*============================================================================*/
bool fis::deFuzzify( void )
{
    bool retValue = false;

    if ( ( nullptr != output ) && ( nullptr != outMF ) ) {
        size_t i;

        for ( i = 0U ; i < nOutputs ; ++i ) {
            deFuzz( &output[ i ] , FIS_DEFUZZ_INIT );
        }
        if ( Mamdani == type ) {
            size_t k;

            for ( k = 0U ; k < nPoints ; ++k ) {
                for ( i = 0U ; i < nOutputs ; ++i ) { /* initialize*/
                    output[ i ].y = 0.0;
                    output[ i ].getNextX( k );
                    output[ i ].value = output[ i ].x;
                }
                fuzzyAggregate();
                for ( i = 0U ; i < nOutputs ; ++i ) {
                    deFuzz( &output[ i ] , FIS_DEFUZZ_COMPUTE );
                    if ( nullptr != output[ i ].xag ) { /*store aggregated*/
                        output[ i ].xag[ k ] = output[ i ].x;
                        output[ i ].yag[ k ] = output[ i ].y;
                    }
                }
            }
        }
        else { /*Sugeno and Tsukamoto systems*/
            for ( i = 0U ; i < nOutputs ; ++i ) { /* initialize*/
                output[ i ].v[ sum_wz ] = 0.0; /*store sum wi*/
                output[ i ].v[ sum_w ] = 0.0; /*store sum zi*wi*/
            }
            fuzzyAggregate();
            for ( i = 0U ; i < nOutputs ; ++i ) { /* initialize*/
                deFuzz( &output[ i ] , FIS_DEFUZZ_COMPUTE );
            }
        }
        for ( i = 0U ; i < nOutputs ; ++i ) {
            output[ i ].value = deFuzz( &output[ i ] , FIS_DEFUZZ_END );
            output[ i ].value = bound( output[ i ].value, output[ i ].min, output[ i ].max );
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool fis::setRuleWeights( real_t *rWeights )
{
    bool retValue = false;
    if ( nullptr != rWeights ) {
        ruleWeight = rWeights;
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool fis::inference( void )
{
    bool retValue = false;

    if ( nullptr != rules  ) {
        size_t i = 0u;

        if ( Q_FIS_RULES_BEGIN == rules[ 0 ] ) {
            inferenceState = &fis::inferenceAntecedent;
            rStrength = 0.0;
            lastConnector = -1;
            ruleCount = 0U;
            i = 1U;
            while ( ( Q_FIS_RULES_END != rules[ i ] ) && ( ruleCount < nRules ) ) {
                i = (this->*inferenceState)( i );
                if ( INFERENCE_ERROR == i ) {
                    break;
                }
                ++i;
            }
        }
        if ( ( Q_FIS_RULES_END == rules[ i ] ) && ( ruleCount == nRules) ) {
            retValue = true;
        }
    }

    return retValue;
}