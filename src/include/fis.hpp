#ifndef QLIBS_FIS
#define QLIBS_FIS

#include "include/types.hpp"

namespace qlibs {

    enum fisShapeMF {
        custommf = 0,   /*!< Custom user-defined Membership function*/
        trimf,          /*!< Triangular Membership function f(a,b,c)*/
        trapmf,         /*!< Trapezoidal Membership function f(a,b,c,d)*/
        gbellmf,        /*!< Generalized bell-shaped Membership function f(a,b,c)*/
        gaussmf,        /*!< Gaussian  Membership function f(s,c)*/
        gauss2mf,       /*!< Gaussian combination Membership function f(s1,c1,s2,c2)*/
        sigmf,          /*!< Sigmoidal  Membership function f(a,c)*/
        dsigmf,         /*!< Difference between two sigmoidal Membership functions f(a1,c1,a2,c2)*/
        psigmf,         /*!< Product of two sigmoidal membership functions f(a1,c1,a2,c2)*/
        pimf,           /*!< Pi-shaped membership function f(a,b,c,d)*/
        smf,            /*!< S-shaped membership function f(a,b)*/
        zmf,            /*!< Z-shaped membership function f(a,b)*/
        singletonmf,    /*!< Singleton Membership Function f(a)*/
        concavemf,      /*!< Concave Membership Function f(i,e)*/
        spikemf,        /*!< Spike Membership Function f(w,c)*/
        linsmf,         /*!< Linear s-shaped saturation membership function f(a,b)*/
        linzmf,         /*!< Linear z-shaped saturation membership function f(a,b)*/
        rectmf,         /*!< Rectangle Membership Function f(s,e)*/
        cosmf,          /*!< Cosine Membership Function f(c,w)*/
        constantmf,     /*!< Constant membership function f(a) [Only for ::Sugeno FIS ]*/
        linearmf,       /*!< Linear membership function f(...) [Only for ::Sugeno FIS ]*/
        tlinsmf,        /*!< Tsukamoto s-shaped saturation membership function f(a,b) [Only for ::Tsukamoto FIS ]*/
        tlinzmf,        /*!< Tsukamoto linzmf membership function f(a,b) [ Only for ::Tsukamoto FIS ]*/
        tconcavemf,     /*!< Tsukamoto z-shaped saturation membership function f(i,e) [Only for ::Tsukamoto FIS ]*/
        tsigmf,         /*!< Tsukamoto Sigmoid membership function f(a,c) [ Only for ::Tsukamoto FIS ]*/
        tsmf,           /*!< Tsukamoto S-Shape membership function f(a,b) [ Only for ::Tsukamoto FIS ]*/
        tzmf,           /*!< Tsukamoto Z-Shape membership function f(a,b) [ Only for ::Tsukamoto FIS ]*/
        /*! @cond  */
        FIS_NUM_MFS        /*!< Number of supported membership functions*/
        /*! @endcond  */
    };

    enum fisDeFuzzMethod {
        centroid = 0,   /*!< Center of gravity of the fuzzy set along the x-axis [ Only for ::Mamdani FIS ]*/
        bisector,       /*!< Vertical line that divides the fuzzy set into two sub-regions of equal area [ Only for ::Mamdani FIS ]**/
        mom,            /*!< Middle of Maximum [ Only for ::Mamdani FIS ]**/
        lom,            /*!< Largest of Maximum [ Only for ::Mamdani FIS ]**/
        som,            /*!< Smallest of Maximum [ Only for ::Mamdani FIS ]**/
        wtaver,         /*!< Weighted average of all rule outputs [ Only for ::Sugeno and ::Tsukamoto FIS ]*/
        wtsum,          /*!< Weighted sum of all rule outputs [ Only for ::Sugeno FIS ]*/
        /*! @cond  */
        FIS_NUM_DEFUZZ      /*!< Number of supported defuzzification methods*/
        /*! @endcond  */
    };

    enum fisParamValue {
        FIS_MIN = 0,   /*!< Minimal value*/
        FIS_PROD,      /*!< Product*/
        FIS_MAX,       /*!< Maximum value*/
        FIS_PROBOR,    /*!< Probabilistic OR*/
        FIS_SUM        /*!< Sum*/
    };

    enum fisParameter {
        FIS_Implication,   /*!< Only ::qFIS_MIN and qFIS_PROD supported*/
        FIS_Aggregation,   /*!< Only ::qFIS_MAX, qFIS_PROBOR and qFIS_SUM supported*/
        FIS_AND,           /*!< Only ::qFIS_MIN and qFIS_PROD supported*/
        FIS_OR,            /*!< Only ::qFIS_MAX and qFIS_PROBOR supported*/
        FIS_EvalPoints     /*!< The number of points for de-fuzzification*/
    };

    enum fisType {
        Mamdani = 0,        /*!< Mamdani inference system. The output of each rule its a fuzzy logic set.*/
        Sugeno,             /*!< Takagi-Sugeno inference system. The output of each rule its a function either linear or constant.*/
        Tsukamoto           /*!< Mamdani inference system. The output of each rule its a fuzzy logic set represented with a monotonic membership function.*/
    };

    enum fisDeFuzzState {
        FIS_DEFUZZ_INIT,
        FIS_DEFUZZ_COMPUTE,
        FIS_DEFUZZ_END
    };

    class fisCore;
    class fis;

    template<fisType fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
    class fisSystem;

    class fisIOBase {
        protected:
            real_t min{ -1.0 };
            real_t max{ 1.0 };
            real_t value{ 0.0 };
        public:
            virtual ~fisIOBase() {};
            fisIOBase() = default;
        friend class fis;
        friend class fisCore;
        template<fisType fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
        friend class fisSystem;
    };

    class fisInput : public fisIOBase {
        public:
            fisInput() = default;
            virtual ~fisInput() {};
        friend class fis;
        friend class fisCore;
        template<fisType fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
        friend class fisSystem;
    };

    class fisOutput : public fisIOBase {
        private:
            fis *owner{ nullptr };
            real_t *xag{ nullptr };
            real_t *yag{ nullptr };
            real_t res{ 0.01 };
            real_t x{ 0.0 };
            real_t y{ 0.0 };
            real_t v[ 4 ] = { 0.0, 0.0, 0.0, 0.0 };
            void getNextX( const size_t i ) noexcept
            {
                x = min + ( ( static_cast<real_t>( i ) + 0.5 )*res );
            }
        public:
            fisOutput() = default;
            virtual ~fisOutput() {};
            bool storeAggregatedRegion( real_t *xData,
                                        real_t *yData,
                                        const size_t n ) noexcept;
            template <size_t numberOfPoints>
            bool storeAggregatedRegion( real_t (&xData)[ numberOfPoints ],
                                        real_t (&yData)[ numberOfPoints ] ) noexcept
            {
                return storeAggregatedRegion( xData, yData, numberOfPoints );
            }
        friend class fis;
        friend class fisCore;
        template<fisType fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
        friend class fisSystem;
    };

    using fisMFFunction = real_t (*)( const fisIOBase * const in,
                                      const real_t *p,
                                      const size_t n );

    class fisCore {
        protected:
            enum defuzzMembers {
                yMax = 0,
                xSmallest = 1,
                xLargest = 2 ,
                sp = 3,
                sum_xy = 0,
                sum_y = 1,
                halfArea = 2,
                currentArea = 3,
                sum_wz = 0,
                sum_w = 1,
            };
            virtual ~fisCore() {}
            fisCore() = default;
            static real_t bound( real_t y,
                                 const real_t minValue,
                                 const real_t maxValue );
            static real_t Min( const real_t a,
                               const real_t b );
            static real_t Prod( const real_t a,
                                const real_t b );
            static real_t Max( const real_t a,
                               const real_t b );
            static real_t ProbOr( const real_t a,
                                  const real_t b );
            static real_t Sum( const real_t a,
                               const real_t b );

            static real_t TriMF( const fisIOBase * const in,
                                 const real_t *p,
                                 const size_t n );
            static real_t TrapMF( const fisIOBase * const in,
                                  const real_t *p,
                                  const size_t n );
            static real_t GBellMF( const fisIOBase * const in,
                                   const real_t *p,
                                   const size_t n );
            static real_t GaussMF( const fisIOBase * const in,
                                   const real_t *p,
                                   const size_t n );
            static real_t Gauss2MF( const fisIOBase * const in,
                                    const real_t *p,
                                    const size_t n );
            static real_t SigMF( const fisIOBase * const in,
                                 const real_t *p,
                                 const size_t n );
            static real_t TSigMF( const fisIOBase * const in,
                                  const real_t *p,
                                  const size_t n );
            static real_t DSigMF( const fisIOBase * const in,
                                  const real_t *p,
                                  const size_t n );
            static real_t PSigMF( const fisIOBase * const in,
                                  const real_t *p,
                                  const size_t n );
            static real_t SMF( const fisIOBase * const in,
                               const real_t *p,
                               const size_t n );
            static real_t TSMF( const fisIOBase * const in,
                                const real_t *p,
                                const size_t n );
            static real_t ZMF( const fisIOBase * const in,
                               const real_t *p,
                               const size_t n );
            static real_t LinSMF( const fisIOBase * const in,
                                  const real_t *p,
                                  const size_t n );
            static real_t LinZMF( const fisIOBase * const in,
                                  const real_t *p,
                                  const size_t n );
            static real_t TZMF( const fisIOBase * const in,
                                const real_t *p,
                                const size_t n );
            static real_t PiMF( const fisIOBase * const in,
                                const real_t *p,
                                const size_t n );
            static real_t SingletonMF( const fisIOBase * const in,
                                       const real_t *p,
                                       const size_t n );
            static real_t ConcaveMF( const fisIOBase * const in,
                                     const real_t *p,
                                     const size_t n );
            static real_t TConcaveMF( const fisIOBase * const in,
                                      const real_t *p,
                                      const size_t n );
            static real_t SpikeMF( const fisIOBase * const in,
                                   const real_t *p,
                                   const size_t n );
            static real_t TLinSMF( const fisIOBase * const in,
                                   const real_t *p,
                                   const size_t n );
            static real_t TLinZMF( const fisIOBase * const in,
                                   const real_t *p,
                                   const size_t n );
            static real_t RectangleMF( const fisIOBase * const in,
                                       const real_t *p,
                                       const size_t n );
            static real_t CosineMF( const fisIOBase * const in,
                                    const real_t *p,
                                    const size_t n );
            static real_t ConstantMF( const fisIOBase * const in,
                                      const real_t *p,
                                      const size_t n );
            static real_t LinearMF( const fisIOBase * const in,
                                    const real_t *p,
                                    const size_t n );

            static real_t deFuzzCentroid( fisOutput * const o,
                                          const fisDeFuzzState stage );
            static real_t deFuzzBisector( fisOutput * const o,
                                          const fisDeFuzzState stage );
            static real_t deFuzzLOM( fisOutput * const o,
                                     const fisDeFuzzState stage );
            static real_t deFuzzSOM( fisOutput * const o,
                                     const fisDeFuzzState stage );
            static real_t deFuzzMOM( fisOutput * const o,
                                     const fisDeFuzzState stage );
            static real_t deFuzzWtAverage( fisOutput * const o,
                                           const fisDeFuzzState stage );
            static real_t deFuzzWtSum( fisOutput * const o,
                                       const fisDeFuzzState stage );
    };

    class fisMF {
        private:
            fisMFFunction shape{ nullptr };
            const real_t *points{ nullptr };
            real_t fx{ 0.0 };
            real_t h{ 0.0 };
            size_t index{ 0U };
            void evalMFAtIndex( fisIOBase *io )
            {
                if ( nullptr != shape ) {
                    fx = h*shape( &io[ index ], points, 1U );
                }
            }
            real_t evalMF( fisIOBase &io )
            {
                real_t y = 0.0;
                if ( nullptr != shape ) {
                    y = h*shape( &io, points, 1U );
                }
                return y;
            }
            real_t evalMF( fisIOBase *io,
                           size_t n )
            {
                real_t y = 0.0;
                if ( nullptr != shape ) {
                    y = shape( io, points, n );
                }
                return y;
            }

        public:
            virtual ~fisMF() {};
            fisMF() = default;
            friend class fis;
    };

    #define FIS_RULE_ITEM_SIZE  1

    #if ( FIS_RULE_ITEM_SIZE == 1 )
        using fisRules = int8_t;
        #define FIS_RULES_MIN_VALUE     INT8_MIN
    #else
        using fisRules = int16_t;
        #define FIS_RULES_MIN_VALUE     INT16_MIN
    #endif

    using fisTag = fisRules;
    /*cstat -MISRAC++2008-0-1-4_b*/
    constexpr fisRules Q_FIS_RULES_BEGIN    = FIS_RULES_MIN_VALUE;
    constexpr fisRules Q_FIS_RULES_END      = FIS_RULES_MIN_VALUE + 1;
    constexpr fisRules Q_FIS_AND            = FIS_RULES_MIN_VALUE + 2;
    constexpr fisRules Q_FIS_OR             = FIS_RULES_MIN_VALUE + 3;
    constexpr fisRules Q_FIS_THEN           = FIS_RULES_MIN_VALUE + 4;
    /*cstat +MISRAC++2008-0-1-4_b*/
    #define Q_FIS_IF_STATEMENT      ,
    #define Q_FIS_AND_STATEMENT     +1),Q_FIS_AND,
    #define Q_FIS_OR_STATEMENT      +1),Q_FIS_OR,
    #define Q_FIS_THEN_STATEMENT    +1),Q_FIS_THEN,
    #define Q_FIS_IS_STATEMENT      ,(
    #define Q_FIS_IS_NOT_STATEMENT  ,-(
    #define Q_FIS_END_STATEMENT     +1)

    #define FIS_RULES_BEGIN         Q_FIS_RULES_BEGIN
    #define FIS_RULES_END           ,Q_FIS_RULES_END
    /** @brief Rule statement to begin a rule sentence */
    #define IF                      Q_FIS_IF_STATEMENT
    /** @brief Rule statement to represent the AND connector */
    #define AND                     Q_FIS_AND_STATEMENT
    /** @brief Rule statement to represent the OR connector */
    #define OR                      Q_FIS_OR_STATEMENT
    /** @brief Rule statement to represent the implication */
    #define THEN                    Q_FIS_THEN_STATEMENT
    /** @brief Rule statement to represent a premise*/
    #define IS                      Q_FIS_IS_STATEMENT
    /** @brief Rule statement to represent a negated premise */
    #define IS_NOT                  Q_FIS_IS_NOT_STATEMENT
    /** @brief Rule statement to end a rule sentence */
    #define END                     Q_FIS_END_STATEMENT


    using fisDeFuzzFunction = real_t (*)( fisOutput * const o, const fisDeFuzzState stage );
    using fuzzyOperator = real_t (*)( const real_t a, const real_t b );

     class fis: public fisCore, private nonCopyable {
        using methods_fcn = real_t (*)( const real_t a, const real_t b );

        template<fisType fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
        friend class fisSystem;
        
        private:
            fisInput *input{ nullptr };
            fisOutput *output{ nullptr };
            fisMF *inMF{ nullptr };
            fisMF *outMF{ nullptr };
            methods_fcn andOp{ &Min };
            methods_fcn orOp{ &Max };
            methods_fcn implicate{ &Min };
            methods_fcn aggregate{ &Max };
            size_t (fis::*inferenceState)( size_t i );
            size_t (fis::*aggregationState)( size_t i );
            fisDeFuzzFunction deFuzz{ &deFuzzCentroid };
            real_t *ruleWeight{ nullptr };
            real_t *wi{ nullptr };
            const fisRules *rules{ nullptr };
            size_t rules_cols{ 0U};
            size_t nInputs{ 0 };
            size_t nOutputs{ 0 };
            size_t nMFInputs{ 0U };
            size_t nMFOutputs{ 0U };
            size_t nPoints{ 100 };
            size_t nRules{ 0U };
            size_t ruleCount{ 0U };
            real_t rStrength{ 0.0 };
            fisRules lastConnector;
            fisType type{ Mamdani };
            bool setMF( fisMF *m,
                        const fisTag io,
                        const fisTag mf,
                        const fisShapeMF s,
                        fisMFFunction customMf,
                        const real_t *cp,
                        const real_t h ) noexcept;
            void evalInputMFs( void ) noexcept;
            void truncateInputs( void ) noexcept;
            static real_t parseFuzzValue( fisMF * const mfIO,
                                          fisRules index ) noexcept;
            fuzzyOperator getFuzzOperator( void ) noexcept;
            size_t inferenceAntecedent( size_t i ) noexcept;
            size_t inferenceReachEnd( size_t i ) noexcept;
            size_t aggregationFindConsequent( size_t i ) noexcept;
            size_t inferenceConsequent( size_t i ) noexcept;
            void fuzzyAggregate( void ) noexcept;
            static const size_t INFERENCE_ERROR;
            fisTag lastTag{ -1 };
        public:
            fis() = default;
            virtual ~fis() {}
            bool setup( const fisType t,
                 fisInput * const inputs,
                 const size_t ni,
                 fisOutput * const outputs,
                 const size_t no,
                 fisMF * const mf_inputs,
                 const size_t nmi,
                 fisMF * const mf_outputs,
                 const size_t nmo,
                 const fisRules * const r,
                 const size_t n,
                 real_t *rWeights = nullptr ) noexcept;
            template <size_t numberInputs, size_t numberOutputs, size_t numberMFinputs, size_t numberMFOutputs, size_t numberRules>
            bool setup( const fisType t,
                 fisInput (&inputs)[ numberInputs ],
                 fisOutput (&outputs)[ numberOutputs ],
                 fisMF (&mf_inputs)[ numberMFinputs ],
                 fisMF (&mf_outputs)[ numberMFOutputs ],
                 const fisRules * const r,
                 real_t (&rWeights)[ numberRules ] ) noexcept
            {
                return setup( t, inputs, numberInputs, outputs, numberOutputs, mf_inputs, numberMFinputs, mf_outputs, numberMFOutputs, r, numberRules, rWeights );
            }
            bool setupInput( const fisTag t,
                             const real_t Min,
                             const real_t Max ) noexcept;
            bool setupOutput( const fisTag t,
                              const real_t Min,
                              const real_t Max ) noexcept;
            bool setInputMF( const fisTag io,
                             const fisTag mf,
                             const fisShapeMF s,
                             const real_t *cp,
                             const real_t h = 1.0 ) noexcept
            {
                return ( ( nullptr != inMF ) && ( nMFInputs > 0U ) ) ? setMF( inMF, io, mf, s, nullptr, cp, h ) : false;
            }
            bool setInputMF( const fisTag io,
                             const fisTag mf,
                             fisMFFunction customMfs,
                             const real_t *cp,
                             const real_t h = 1.0 ) noexcept
            {
                return ( ( nullptr != inMF ) && ( nMFInputs > 0U ) ) ? setMF( inMF, io, mf, custommf, customMfs, cp, h ) : false;
            }
            bool setOutputMF( const fisTag io,
                              const fisTag mf,
                              const fisShapeMF s,
                              const real_t *cp,
                              const real_t h = 1.0 ) noexcept
            {
                return ( ( nullptr != outMF ) && ( nMFOutputs > 0U ) ) ? setMF( outMF, io, mf, s, nullptr, cp, h ) : false; 
            }
            bool setOutputMF( const fisTag io,
                              const fisTag mf,
                              fisMFFunction customMfs,
                              const real_t *cp,
                              const real_t h = 1.0 ) noexcept
            {
                return ( ( nullptr != outMF ) && ( nMFOutputs > 0U ) ) ? setMF( outMF, io, mf, custommf, customMfs, cp, h ) : false; 
            }
            bool setInput( const fisTag t,
                           const real_t value ) noexcept;
            bool getOutput( const fisTag t,
                            real_t &value ) const noexcept;
            bool setParameter( const fisParameter p,
                               const fisParamValue x ) noexcept;
            bool setDeFuzzMethod( fisDeFuzzMethod m ) noexcept;
            bool fuzzify( void ) noexcept;
            bool deFuzzify( void ) noexcept;
            bool inference( void ) noexcept;
            bool setRuleWeights( real_t *rWeights ) noexcept;
            size_t getNumberOfPoints( void ) const noexcept
            {
                return nPoints;
            }
            real_t operator[]( fisTag outTag ) const
            {
                /*cstat -CERT-STR34-C*/
                return ( static_cast<size_t>( outTag ) < nOutputs ) ? output[ outTag ].value :  input[ nOutputs -1 ].value;
                /*cstat +CERT-STR34-C*/
            }

            fis& operator<<( const fisTag& tag ) {
                if ( static_cast<size_t>( tag ) < nInputs ) {
                    lastTag = tag;
                }
                return *this;
            }
            fis& operator<<( const int& value ) {
                if ( lastTag >= 0 ) {
                    input[ lastTag ].value = static_cast<real_t>(value);
                }
                return *this;
            }
            fis& operator<<( const real_t& value ) {
                if ( lastTag >= 0 ) {
                    input[ lastTag ].value = value;
                }
                return *this;
            }

        friend class fisCore;
    };

    template<fisType fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
    class fisSystem {
        private:
            fis sys;
            const fisRules *rules;
            fisInput inputs[ numberOfInputs ];
            fisOutput outputs[ numberOfOutputs ];
            fisMF MFin[ numberOfInputSets ], MFout[ numberOfOutputSets ];
            real_t ruleStrength[ numberOfRules ];
        public:
            constexpr fisSystem( const fisRules *Rules ) : rules( Rules ) {};
            inline bool setup( void )
            {
                return sys.setup( fType,
                                  inputs, numberOfInputs,
                                  outputs, numberOfOutputs,
                                  MFin, numberOfInputSets,
                                  MFout, numberOfOutputSets,
                                  rules, numberOfRules,
                                  ruleStrength );
            }
            inline bool setupInput( const fisTag t,
                                    const real_t Min,
                                    const real_t Max ) noexcept
            {
                return sys.setupInput( t, Min, Max );
            }
            inline bool setupOutput( const fisTag t,
                                     const real_t Min,
                                     const real_t Max ) noexcept
            {
                return sys.setupOutput( t, Min, Max );
            }
            inline bool setInputMF( const fisTag io,
                                    const fisTag mf,
                                    const fisShapeMF s,
                                    const real_t *cp,
                                    const real_t h = 1.0 ) noexcept
            {
                return sys.setInputMF( io, mf, s, cp, h );
            }
            inline bool setInputMF( const fisTag io,
                                    const fisTag mf,
                                    fisMFFunction customMfs,
                                    const real_t *cp,
                                    const real_t h = 1.0 ) noexcept
            {
                return sys.setInputMF( io, mf, customMfs, cp, h );
            }
            inline bool setOutputMF( const fisTag io,
                                     const fisTag mf,
                                     const fisShapeMF s,
                                     const real_t *cp,
                                     const real_t h = 1.0 ) noexcept
            {
                return sys.setOutputMF( io, mf, s, cp, h );
            }
            inline bool setOutputMF( const fisTag io,
                                     const fisTag mf,
                                     fisMFFunction customMfs,
                                     const real_t *cp,
                                     const real_t h = 1.0 ) noexcept
            {
                return sys.setOutputMF( io, mf, customMfs, cp, h );
            }
            inline bool fuzzify( void ) noexcept
            {
                return sys.fuzzify();
            }
            inline bool deFuzzify( void ) noexcept
            {
                return sys.deFuzzify();
            }
            inline bool inference( void ) noexcept
            {
                return sys.inference();
            }

            inline bool setInput( const fisTag t,
                                  const real_t value ) noexcept
            {
                return sys.setInput( t, value );
            }
            inline bool getOutput( const fisTag t,
                                   real_t &value ) const noexcept
            {
                return sys.getOutput( t, value );
            }
            inline bool setParameter( const fisParameter p,
                                      const fisParamValue x ) noexcept
            {
                return sys.setParameter( p, x );
            }
            inline bool setDeFuzzMethod( fisDeFuzzMethod m ) noexcept
            {
                return sys.setDeFuzzMethod( m );
            }
            inline bool setRuleWeights( real_t *rWeights ) noexcept
            {
                return sys.setRuleWeights( rWeights );
            }
            size_t getNumberOfPoints( void ) const noexcept
            {
                return sys.nPoints;
            }

            real_t operator[]( fisTag outTag ) const
            {
                /*cstat -CERT-STR34-C*/
                return ( static_cast<size_t>( outTag ) < sys.nOutputs ) ? sys.output[ outTag ].value :  sys.input[ sys.nOutputs -1 ].value;
                /*cstat +CERT-STR34-C*/
            }

            fisSystem& operator<<( const fisTag& tag ) {
                if ( static_cast<size_t>( tag ) < sys.nInputs ) {
                    sys.lastTag = tag;
                }
                return *this;
            }
            fisSystem& operator<<( const int& value ) {
                if ( sys.lastTag >= 0 ) {
                    sys.input[ sys.lastTag ].value = static_cast<real_t>( value );
                }
                return *this;
            }
            fisSystem& operator<<( const real_t& value ) {
                if ( sys.lastTag >= 0 ) {
                    sys.input[ sys.lastTag ].value = value;
                }
                return *this;
            }
    };


}

#endif /*QLIBS_FIS*/