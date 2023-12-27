/*!
 * @file fis.hpp
 * @author J. Camilo Gomez C.
 * @version 1.0.1
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Fuzzy Inference System (FIS) Engine
 **/

#ifndef QLIBS_FIS
#define QLIBS_FIS

#include "include/types.hpp"

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /** @addtogroup qfis FIS - Fuzzy Inference System
    * @brief Fuzzy Inference System (FIS) Engine
    *  @{
    */

    /**
    * @brief An enum with all the possible values to specify a membership
    * function.
    */
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

    /**
    * @brief An enum with all the possible de-Fuzzyfication methods.
    */
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

    /**
    * @brief An enum with the supported parameter values
    */
    enum fisParamValue {
        FIS_MIN = 0,   /*!< Minimal value*/
        FIS_PROD,      /*!< Product*/
        FIS_MAX,       /*!< Maximum value*/
        FIS_PROBOR,    /*!< Probabilistic OR*/
        FIS_SUM        /*!< Sum*/
    };

    /**
    * @brief An enum with the allowed parameters that can be set on a FIS instance
    */
    enum fisParameter {
        FIS_Implication,   /*!< Only ::FIS_MIN and FIS_PROD supported*/
        FIS_Aggregation,   /*!< Only ::FIS_MAX, FIS_PROBOR and qFIS_SUM supported*/
        FIS_AND,           /*!< Only ::FIS_MIN and FIS_PROD supported*/
        FIS_OR,            /*!< Only ::FIS_MAX and FIS_PROBOR supported*/
        FIS_EvalPoints     /*!< The number of points for de-fuzzification*/
    };

    /**
    * @brief An enum with the inference system types supported by qFIS
    */
    enum fisType {
        Mamdani = 0,        /*!< Mamdani inference system. The output of each rule its a fuzzy logic set.*/
        Sugeno,             /*!< Takagi-Sugeno inference system. The output of each rule its a function either linear or constant.*/
        Tsukamoto           /*!< Mamdani inference system. The output of each rule its a fuzzy logic set represented with a monotonic membership function.*/
    };

    /*! @cond  */
    enum fisDeFuzzState {
        FIS_DEFUZZ_INIT,
        FIS_DEFUZZ_COMPUTE,
        FIS_DEFUZZ_END
    };
    /*! @endcond  */
    
    /*! @cond  */
    class fisCore;
    class fis;
    template<fisType fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
    class fisSystem;
    /*! @endcond  */

    /*! @cond  */
    class fisIOBase {
        protected:
            /*! @cond  */
            real_t min{ -1.0 };
            real_t max{ 1.0 };
            real_t value{ 0.0 };
            /*! @endcond  */
        public:
            virtual ~fisIOBase() {};
            fisIOBase() = default;
        friend class fis;
        friend class fisCore;
        template<fisType fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
        friend class fisSystem;
    };
    /*! @endcond  */

    /**
    * @brief A FIS Input object
    * @details The instance should be initialized using the fis::inputSetup() method.
    */
    class fisInput : public fisIOBase {
        public:
            fisInput() = default;
            virtual ~fisInput() {};
        friend class fis;
        friend class fisCore;
        template<fisType fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
        friend class fisSystem;
    };

    /**
    * @brief A FIS Output object
    * @details The instance should be initialized using the fis::outputSetup() method.
    */
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

    /*! @cond  */
    using fisMFFunction = real_t (*)( const fisIOBase * const in,
                                      const real_t *p,
                                      const size_t n );
    /*! @endcond  */

    class fisCore {
        protected:
            /*! @cond  */
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
            /*! @endcond  */
    };

    /**
    * @brief A FIS Membership Function
    * @details The instance should be initialized using the fis::setInputMF() 
    * or fis::setOutputMF() methods.
    */
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
        /**
        * @brief Type definition to instantiate a set of fuzzy rules
        * @details Rules are defined by combining I/O and membership function tags
        * with the following statements:
        * 
        * #FIS_RULES_BEGIN, #IF, #IS, #IS_NOT, #AND, #OR, #THEN, #END and
        * #FIS_RULES_END
        * 
        * Example:
        * @code{.c}
        * static const fisRules rules[] = { 
        *     FIS_RULES_BEGIN
        *         IF service IS service_poor OR food IS food_rancid THEN tip IS tip_cheap END
        *         IF service IS service_good THEN tip IS tip_average END
        *         IF service IS service_excellent OR food IS food_delicious THEN tip IS tip_generous END
        *     FIS_RULES_END
        * };
        * @endcode
        */
        /*! @cond  */
        using fisRules = int8_t;
        #define FIS_RULES_MIN_VALUE     INT8_MIN
        /*! @endcond  */
    #else
        /**
        * @brief Type definition to instantiate a set of fuzzy rules
        * @details Rules are defined by combining I/O and membership function tags
        * with the following statements:
        * 
        * #FIS_RULES_BEGIN, #IF, #IS, #IS_NOT, #AND, #OR, #THEN, #END and
        * #FIS_RULES_END
        * 
        * Example:
        * @code{.c}
        * static const fisRules rules[] = { 
        *     FIS_RULES_BEGIN
        *         IF service IS service_poor OR food IS food_rancid THEN tip IS tip_cheap END
        *         IF service IS service_good THEN tip IS tip_average END
        *         IF service IS service_excellent OR food IS food_delicious THEN tip IS tip_generous END
        *     FIS_RULES_END
        * };
        * @endcode
        */
        /*! @cond  */
        using fisRules = int16_t;
        #define FIS_RULES_MIN_VALUE     INT16_MIN
        /*! @endcond  */
    #endif

    using fisTag = fisRules;
    /*cstat -MISRAC++2008-0-1-4_b*/
    /*! @cond  */
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
    /*! @endcond  */

    /**
    * @brief Start a Fuzzy rule set.
    * The #FIS_RULES_BEGIN statement is used to declare the starting point of
    * a FIS rule set. It should be placed at the start of the rules enumeration.
    * #FIS_RULES_END declare the end of the FIS rule set.
    * @see #FIS_RULES_END
    * @warning Only one segment is allowed inside a fuzzy rule set.
    * @note It must always be used together with a matching #FIS_RULES_END 
    * statement.
    * Example:
    * @code{.c}
    * static const fisRules rules[] = { 
    *     FIS_RULES_BEGIN
    *       
    *     FIS_RULES_END
    * };
    * @endcode
    */
    #define FIS_RULES_BEGIN         Q_FIS_RULES_BEGIN
    /**
    * @brief Ends a Fuzzy rule set.
    * The #FIS_RULES_END statement is used to finalize the declaration of a 
    * FIS rule set. It should be placed at the end of the rules enumeration.
    * #FIS_RULES_BEGIN declare the start point of the FIS rule set.
    * @see #FIS_RULES_BEGIN
    * @warning Only one segment is allowed inside a fuzzy rule set.
    * @note It must always be used together with a matching #FIS_RULES_BEGIN 
    * statement.
    * Example:
    * @code{.c}
    * static const fisRules_t rules[] = { 
    *     FIS_RULES_BEGIN
    *       
    *     FIS_RULES_END
    * };
    * @endcode
    */
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

    /**
    * @brief A FIS(Fuzzy Inference System) object
    * @details The instance should be initialized using the fis::setup() method.
    */
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

            /**
            * @brief Setup and initialize the FIS instance.
            * @note Default configuration : AND = Min, OR = Max, Implication = Min
            * Aggregation = Max, EvalPoints = 100
            * @param[in] t Type of inference ::Mamdani, ::Sugeno or ::Tsukamoto.
            * @param[in] inputs An array with all the system inputs as fisInput
            * objects.
            * @param[in] ni The number of elements in the @a inputs array.
            * @param[in] outputs An array with all the system outputs as fisOutput
            * objects.
            * @param[in] no The number of elements in the @a outputs array.
            * @param[in] mf_inputs An array with all the membership functions related to
            * the inputs. This should be an array of fisMF objects.
            * @param[in] nmi The number of elements in the @a mf_inputs array.
            * @param[in] mf_outputs An array with all the membership functions related to
            * the outputs. This should be an array of fisMF objects.
            * @param[in] nmo The number of elements in the @a mf_outputs array.
            * @param[in] r The rules set.
            * @param[in] n Number of rules
            * @param[in] rWeights An array of size @a n were the rule strengths will be stored.
            * @return @c true on success, otherwise return @c false.
            */
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

            /**
            * @brief Setup and initialize the FIS instance.
            * @note Default configuration : AND = Min, OR = Max, Implication = Min
            * Aggregation = Max, EvalPoints = 100
            * @param[in] t Type of inference ::Mamdani, ::Sugeno or ::Tsukamoto.
            * @param[in] inputs An array with all the system inputs as fisInput
            * objects.
            * @param[in] outputs An array with all the system outputs as fisOutput
            * objects.
            * @param[in] mf_inputs An array with all the membership functions related to
            * the inputs. This should be an array of fisMF objects.
            * @param[in] mf_outputs An array with all the membership functions related to
            * the outputs. This should be an array of fisMF objects.
            * @param[in] r The rules set.
            * @param[in] rWeights An array of size @a n were the rule strengths will be stored.
            * @return @c true on success, otherwise return @c false.
            */
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

            /**
            * @brief Setup the input with the specified tag and set limits for it
            * @param[in] t The input tag
            * @param[in] min Minimum allowed value for this input
            * @param[in] max Max allowed value for this input
            * @return @c true on success, otherwise return @c false.
            */
            bool setupInput( const fisTag t,
                             const real_t Min,
                             const real_t Max ) noexcept;
            /**
            * @brief Setup the output with the specified tag and set limits for it
            * @param[in] t The output tag
            * @param[in] min Minimum allowed value for this output
            * @param[in] max Max allowed value for this output
            * @return @c true on success, otherwise return @c false.
            */
            bool setupOutput( const fisTag t,
                              const real_t Min,
                              const real_t Max ) noexcept;

            /**
            * @brief Set the input tag and points for the specified membership
            * function
            * @param[in] io The input tag related with this membership function
            * @param[in] mf The user-defined tag for this membership function
            * @param[in] s The wanted shape/form for this membership function, can
            * be one of the following: ::trimf, ::trapmf, ::gbellmf, ::gaussmf,
            * ::gauss2mf, ::sigmf, ::dsigmf, ::psigmf, ::pimf, ::smf, ::zmf,
            * ::singletonmf, ::concavemf, ::spikemf, ::linsmf, ::linzmf, ::rectmf,
            * ::cosmf.
            * @note For ::Sugeno FIS, an output membership function should be one of the
            * following: ::constantmf, ::linearmf.
            * @note For ::Tsukamoto FIS, an output membership function should be one the
            * following monotonic functions : ::tlinsmf, ::tlinzmf, ::tsmf, ::tzmf, 
            * ::tconcavemf
            * @param[in] cp Points or coefficients of the membership function.
            * @param[in] h Height of the membership function.
            * @note Height parameter @a h does not apply for output membership functions
            * on ::Sugeno and ::Tsukamoto inference systems. [ 0 <= h <= 1]
            * @return @c true on success, otherwise return @c false.
            */
            bool setInputMF( const fisTag io,
                             const fisTag mf,
                             const fisShapeMF s,
                             const real_t *cp,
                             const real_t h = 1.0 ) noexcept
            {
                return ( ( nullptr != inMF ) && ( nMFInputs > 0U ) ) ? setMF( inMF, io, mf, s, nullptr, cp, h ) : false;
            }

            /**
            * @brief Set the input tag and points for the specified membership
            * function
            * @param[in] io The input tag related with this membership function
            * @param[in] mf The user-defined tag for this membership function
            * @param[in] customMfs Custom user-defined membership function.
            * @param[in] cp Points or coefficients of the membership function.
            * @param[in] h Height of the membership function.
            * @note Height parameter @a h does not apply for output membership functions
            * on ::Sugeno and ::Tsukamoto inference systems. [ 0 <= h <= 1]
            * @return @c true on success, otherwise return @c false.
            */
            bool setInputMF( const fisTag io,
                             const fisTag mf,
                             fisMFFunction customMfs,
                             const real_t *cp,
                             const real_t h = 1.0 ) noexcept
            {
                return ( ( nullptr != inMF ) && ( nMFInputs > 0U ) ) ? setMF( inMF, io, mf, custommf, customMfs, cp, h ) : false;
            }

            /**
            * @brief Set the output tag and points for the specified membership
            * function
            * @param[in] io The output tag related with this membership function
            * @param[in] mf The user-defined tag for this membership function
            * @param[in] s The wanted shape/form for this membership function, can
            * be one of the following: ::trimf, ::trapmf, ::gbellmf, ::gaussmf,
            * ::gauss2mf, ::sigmf, ::dsigmf, ::psigmf, ::pimf, ::smf, ::zmf,
            * ::singletonmf, ::concavemf, ::spikemf, ::linsmf, ::linzmf, ::rectmf,
            * ::cosmf.
            * @note For ::Sugeno FIS, an output membership function should be one of the
            * following: ::constantmf, ::linearmf.
            * @note For ::Tsukamoto FIS, an output membership function should be one the
            * following monotonic functions : ::tlinsmf, ::tlinzmf, ::tsmf, ::tzmf, 
            * ::tconcavemf
            * @param[in] cp Points or coefficients of the membership function.
            * @param[in] h Height of the membership function.
            * @note Height parameter @a h does not apply for output membership functions
            * on ::Sugeno and ::Tsukamoto inference systems. [ 0 <= h <= 1]
            * @return @c true on success, otherwise return @c false.
            */
            bool setOutputMF( const fisTag io,
                              const fisTag mf,
                              const fisShapeMF s,
                              const real_t *cp,
                              const real_t h = 1.0 ) noexcept
            {
                return ( ( nullptr != outMF ) && ( nMFOutputs > 0U ) ) ? setMF( outMF, io, mf, s, nullptr, cp, h ) : false; 
            }

            /**
            * @brief Set the output tag and points for the specified membership
            * function
            * @param[in] io The output tag related with this membership function
            * @param[in] mf The user-defined tag for this membership function
            * @param[in] customMfs Custom user-defined membership function.
            * @param[in] cp Points or coefficients of the membership function.
            * @param[in] h Height of the membership function.
            * @note Height parameter @a h does not apply for output membership functions
            * on ::Sugeno and ::Tsukamoto inference systems. [ 0 <= h <= 1]
            * @return @c true on success, otherwise return @c false.
            */
            bool setOutputMF( const fisTag io,
                              const fisTag mf,
                              fisMFFunction customMfs,
                              const real_t *cp,
                              const real_t h = 1.0 ) noexcept
            {
                return ( ( nullptr != outMF ) && ( nMFOutputs > 0U ) ) ? setMF( outMF, io, mf, custommf, customMfs, cp, h ) : false; 
            }

            /**
            * @brief Set a crisp value of the input with the specified tag.
            * @param[in] t The input tag
            * @param[in] value The crisp value to set
            * @return @c true on success, otherwise return @c false.
            */
            bool setInput( const fisTag t,
                           const real_t value ) noexcept;

            /**
            * @brief Get the de-fuzzified crisp value from the output with the
            * specified tag.
            * @param[in] t The output tag
            * @param[out] value The variable where the output will be stored
            * @return The requested de-fuzzified crisp value.
            */
            bool getOutput( const fisTag t,
                            real_t &value ) const noexcept;

            /**
            * @brief Set parameters of the FIS instance.
            * @param[in] p The requested parameter to change/set.
            * @param[in] x The value of the parameter to set.
            * @return @c true on success, otherwise return @c false.
            */
            bool setParameter( const fisParameter p,
                               const fisParamValue x ) noexcept;

            /**
            * @brief Change the default de-Fuzzification method of the FIS instance.
            * @param[in] m The de-fuzzification method: use one of the following :
            *  ::centroid, ::bisector, ::mom, ::lom, ::som, ::wtaver, ::wtsum
            * @note ::centroid, ::bisector, ::mom, ::lom and ::som only apply for a
            * ::Mamdani FIS
            * @note ::wtaver and ::wtsum only apply for a ::Sugeno FIS.
            * @note ::wtaver only apply for a ::Tsukamoto FIS
            * @return @c true on success, otherwise return @c false
            */
            bool setDeFuzzMethod( fisDeFuzzMethod m ) noexcept;

            /**
            * @brief Perform the fuzzification operation over the crisp inputs on the
            * requested FIS object
            * @pre I/Os and fuzzy sets must be previously initialized by fis::setupInput(),
            * fis::setupOutput(), fis::setInputMF(), fis::setOutputMF and fis::setup() respectively.
            * @return @c true on success, otherwise return @c false.
            */
            bool fuzzify( void ) noexcept;

            /**
            * @brief Perform the de-Fuzzification operation to compute the crisp outputs.
            * @pre The instance should have already invoked the inference process 
            * successfully with fis::inference()
            * @note By default, this method, uses the Centroid method on
            * ::Mamdani type FIS and weight-average on ::Sugeno type FIS. To change
            * the default settings use the fis::setDeFuzzMethod() function.
            * @return @c true on success, otherwise return @c false.
            */
            bool deFuzzify( void ) noexcept;

            /**
            * @brief Perform the inference process on the FIS object
            * @pre The instance should have already invoked the fuzzification operation 
            * successfully with fis::fuzzify()
            * @return @c true on success, otherwise return @c false.
            */
            bool inference( void ) noexcept;

            /**
            * @brief Set weights to the rules of the inference system.
            * @pre I/Os and fuzzy sets must be previously initialized by fis::setupInput(),
            * fis::setupOutput(), fis::setInputMF(), fis::setOutputMF and fis::setup() respectively.
            * @param[in] rWeights An array with the values of every rule weight;
            * @return @c true on success, otherwise return @c false.
            */
            bool setRuleWeights( real_t *rWeights ) noexcept;

            /**
            * @brief Get the number of points used on Mamdani to perform the 
            * de-fuzzification proccess
            * @return The number for points for de-fuzzification.
            */
            size_t getNumberOfPoints( void ) const noexcept
            {
                return nPoints;
            }

            /**
            * @brief Get the de-fuzzified crisp value from the output with the
            * specified tag.
            * @param[in] outTag The output tag
            * @return The requested de-fuzzified crisp value.
            */
            real_t operator[]( fisTag outTag ) const
            {
                /*cstat -CERT-STR34-C*/
                return ( static_cast<size_t>( outTag ) < nOutputs ) ? output[ outTag ].value :  input[ nOutputs -1 ].value;
                /*cstat +CERT-STR34-C*/
            }

            /**
            * @brief Select the input to set using with the specified tag.
            * @param[in] tag The input tag
            * @return The fis object you invoked << upon.
            */
            fis& operator<<( const fisTag& tag ) {
                if ( static_cast<size_t>( tag ) < nInputs ) {
                    lastTag = tag;
                }
                return *this;
            }

            /**
            * @brief The value to set the previously selected input.
            * @param[in] value The crisp value to set
            * @return The fis object you invoked << upon.
            */
            fis& operator<<( const int& value ) {
                if ( lastTag >= 0 ) {
                    input[ lastTag ].value = static_cast<real_t>(value);
                }
                return *this;
            }

            /**
            * @brief The value to set the previously selected input.
            * @param[in] value The crisp value to set
            * @return The fis object you invoked << upon.
            */
            fis& operator<<( const real_t& value ) {
                if ( lastTag >= 0 ) {
                    input[ lastTag ].value = value;
                }
                return *this;
            }

        friend class fisCore;
    };

    /**
    * @brief A wrapper for the FIS object
    * @details The instance should be initialized using the fisSystem::setup() method.
    * @tparam fType Type of inference ::Mamdani, ::Sugeno or ::Tsukamoto.
    * @tparam numberOfInputs The number of inputs of the FIS system.
    * @tparam numberOfOutputs The number of outputs of the FIS system.
    * @tparam numberOfInputSets The number of sets/membership functions for 
    * the inputs.
    * @tparam numberOfOutputSets The number of sets/membership functions for 
    * the outputs.
    * @tparam numberOfRules Number of rules
    */
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

            /**
            * @brief Setup and initialize the FIS instance.
            * @return @c true on success, otherwise return @c false.
            */
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

            /**
            * @brief Setup the input with the specified tag and set limits for it
            * @param[in] t The input tag
            * @param[in] Min Minimum allowed value for this input
            * @param[in] Max Max allowed value for this input
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setupInput( const fisTag t,
                                    const real_t Min,
                                    const real_t Max ) noexcept
            {
                return sys.setupInput( t, Min, Max );
            }

            /**
            * @brief Setup the output with the specified tag and set limits for it
            * @param[in] t The output tag
            * @param[in] Min Minimum allowed value for this output
            * @param[in] Max Max allowed value for this output
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setupOutput( const fisTag t,
                                     const real_t Min,
                                     const real_t Max ) noexcept
            {
                return sys.setupOutput( t, Min, Max );
            }

            /**
            * @brief Set the input tag and points for the specified membership
            * function
            * @param[in] io The input tag related with this membership function
            * @param[in] mf The user-defined tag for this membership function
            * @param[in] s The wanted shape/form for this membership function, can
            * be one of the following: ::trimf, ::trapmf, ::gbellmf, ::gaussmf,
            * ::gauss2mf, ::sigmf, ::dsigmf, ::psigmf, ::pimf, ::smf, ::zmf,
            * ::singletonmf, ::concavemf, ::spikemf, ::linsmf, ::linzmf, ::rectmf,
            * ::cosmf.
            * @note For ::Sugeno FIS, an output membership function should be one of the
            * following: ::constantmf, ::linearmf.
            * @note For ::Tsukamoto FIS, an output membership function should be one the
            * following monotonic functions : ::tlinsmf, ::tlinzmf, ::tsmf, ::tzmf, 
            * ::tconcavemf
            * @param[in] cp Points or coefficients of the membership function.
            * @param[in] h Height of the membership function.
            * @note Height parameter @a h does not apply for output membership functions
            * on ::Sugeno and ::Tsukamoto inference systems. [ 0 <= h <= 1]
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setInputMF( const fisTag io,
                                    const fisTag mf,
                                    const fisShapeMF s,
                                    const real_t *cp,
                                    const real_t h = 1.0 ) noexcept
            {
                return sys.setInputMF( io, mf, s, cp, h );
            }

            /**
            * @brief Set the input tag and points for the specified membership
            * function
            * @param[in] io The input tag related with this membership function
            * @param[in] mf The user-defined tag for this membership function
            * @param[in] customMfs Custom user-defined membership function.
            * @param[in] cp Points or coefficients of the membership function.
            * @param[in] h Height of the membership function.
            * @note Height parameter @a h does not apply for output membership functions
            * on ::Sugeno and ::Tsukamoto inference systems. [ 0 <= h <= 1]
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setInputMF( const fisTag io,
                                    const fisTag mf,
                                    fisMFFunction customMfs,
                                    const real_t *cp,
                                    const real_t h = 1.0 ) noexcept
            {
                return sys.setInputMF( io, mf, customMfs, cp, h );
            }

            /**
            * @brief Set the output tag and points for the specified membership
            * function
            * @param[in] io The output tag related with this membership function
            * @param[in] mf The user-defined tag for this membership function
            * @param[in] s The wanted shape/form for this membership function, can
            * be one of the following: ::trimf, ::trapmf, ::gbellmf, ::gaussmf,
            * ::gauss2mf, ::sigmf, ::dsigmf, ::psigmf, ::pimf, ::smf, ::zmf,
            * ::singletonmf, ::concavemf, ::spikemf, ::linsmf, ::linzmf, ::rectmf,
            * ::cosmf.
            * @note For ::Sugeno FIS, an output membership function should be one of the
            * following: ::constantmf, ::linearmf.
            * @note For ::Tsukamoto FIS, an output membership function should be one the
            * following monotonic functions : ::tlinsmf, ::tlinzmf, ::tsmf, ::tzmf, 
            * ::tconcavemf
            * @param[in] cp Points or coefficients of the membership function.
            * @param[in] h Height of the membership function.
            * @note Height parameter @a h does not apply for output membership functions
            * on ::Sugeno and ::Tsukamoto inference systems. [ 0 <= h <= 1]
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setOutputMF( const fisTag io,
                                     const fisTag mf,
                                     const fisShapeMF s,
                                     const real_t *cp,
                                     const real_t h = 1.0 ) noexcept
            {
                return sys.setOutputMF( io, mf, s, cp, h );
            }

            /**
            * @brief Set the output tag and points for the specified membership
            * function
            * @param[in] io The output tag related with this membership function
            * @param[in] mf The user-defined tag for this membership function
            * @param[in] customMfs Custom user-defined membership function.
            * @param[in] cp Points or coefficients of the membership function.
            * @param[in] h Height of the membership function.
            * @note Height parameter @a h does not apply for output membership functions
            * on ::Sugeno and ::Tsukamoto inference systems. [ 0 <= h <= 1]
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setOutputMF( const fisTag io,
                                     const fisTag mf,
                                     fisMFFunction customMfs,
                                     const real_t *cp,
                                     const real_t h = 1.0 ) noexcept
            {
                return sys.setOutputMF( io, mf, customMfs, cp, h );
            }

            /**
            * @brief Perform the fuzzification operation over the crisp inputs on the
            * requested FIS object
            * @pre I/Os and fuzzy sets must be previously initialized by fis::setupInput(),
            * fisSystem::setupOutput(), fisSystem::setInputMF(), fisSystem::setOutputMF 
            * and fisSystem::setup() respectively.
            * @return @c true on success, otherwise return @c false.
            */
            inline bool fuzzify( void ) noexcept
            {
                return sys.fuzzify();
            }

            /**
            * @brief Perform the de-Fuzzification operation to compute the crisp outputs.
            * @pre The instance should have already invoked the inference process 
            * successfully with fisSystem::inference()
            * @note By default, this method, uses the Centroid method on
            * ::Mamdani type FIS and weight-average on ::Sugeno type FIS. To change
            * the default settings use the fis::setDeFuzzMethod() function.
            * @return @c true on success, otherwise return @c false.
            */
            inline bool deFuzzify( void ) noexcept
            {
                return sys.deFuzzify();
            }

            /**
            * @brief Perform the inference process on the FIS object
            * @pre The instance should have already invoked the fuzzification operation 
            * successfully with fisSystem::fuzzify()
            * @return @c true on success, otherwise return @c false.
            */
            inline bool inference( void ) noexcept
            {
                return sys.inference();
            }

            /**
            * @brief Set a crisp value of the input with the specified tag.
            * @param[in] t The input tag
            * @param[in] value The crisp value to set
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setInput( const fisTag t,
                                  const real_t value ) noexcept
            {
                return sys.setInput( t, value );
            }

            /**
            * @brief Get the de-fuzzified crisp value from the output with the
            * specified tag.
            * @param[in] t The output tag
            * @param[out] value The variable where the output will be stored.
            * @return The requested de-fuzzified crisp value.
            */
            inline bool getOutput( const fisTag t,
                                   real_t &value ) const noexcept
            {
                return sys.getOutput( t, value );
            }

            /**
            * @brief Set parameters of the FIS instance.
            * @param[in] p The requested parameter to change/set.
            * @param[in] x The value of the parameter to set.
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setParameter( const fisParameter p,
                                      const fisParamValue x ) noexcept
            {
                return sys.setParameter( p, x );
            }

            /**
            * @brief Change the default de-Fuzzification method of the FIS instance.
            * @param[in] m The de-fuzzification method: use one of the following :
            *  ::centroid, ::bisector, ::mom, ::lom, ::som, ::wtaver, ::wtsum
            * @note ::centroid, ::bisector, ::mom, ::lom and ::som only apply for a
            * ::Mamdani FIS
            * @note ::wtaver and ::wtsum only apply for a ::Sugeno FIS.
            * @note ::wtaver only apply for a ::Tsukamoto FIS
            * @return @c true on success, otherwise return @c false
            */
            inline bool setDeFuzzMethod( fisDeFuzzMethod m ) noexcept
            {
                return sys.setDeFuzzMethod( m );
            }

            /**
            * @brief Set weights to the rules of the inference system.
            * @pre I/Os and fuzzy sets must be previously initialized by 
            * fisSystem::setupInput(), fisSystem::setupOutput(), 
            * fisSystem::setInputMF(), fisSystem::setOutputMF and 
            * fisSystem::setup() respectively.
            * @param[in] rWeights An array with the values of every rule weight;
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setRuleWeights( real_t *rWeights ) noexcept
            {
                return sys.setRuleWeights( rWeights );
            }

            /**
            * @brief Get the number of points used on Mamdani to perform the 
            * de-fuzzification proccess
            * @return The number for points for de-fuzzification.
            */
            size_t getNumberOfPoints( void ) const noexcept
            {
                return sys.nPoints;
            }

            /**
            * @brief Get the de-fuzzified crisp value from the output with the
            * specified tag.
            * @param[in] outTag The output tag
            * @return The requested de-fuzzified crisp value.
            */
            real_t operator[]( fisTag outTag ) const
            {
                /*cstat -CERT-STR34-C*/
                return ( static_cast<size_t>( outTag ) < sys.nOutputs ) ? sys.output[ outTag ].value :  sys.input[ sys.nOutputs -1 ].value;
                /*cstat +CERT-STR34-C*/
            }

            /**
            * @brief Select the input to set using with the specified tag.
            * @param[in] tag The input tag
            * @return The fisSystem object you invoked << upon.
            */
            fisSystem& operator<<( const fisTag& tag ) {
                if ( static_cast<size_t>( tag ) < sys.nInputs ) {
                    sys.lastTag = tag;
                }
                return *this;
            }

            /**
            * @brief The value to set the previously selected input.
            * @param[in] value The crisp value to set
            * @return The fisSystem object you invoked << upon.
            */
            fisSystem& operator<<( const int& value ) {
                if ( sys.lastTag >= 0 ) {
                    sys.input[ sys.lastTag ].value = static_cast<real_t>( value );
                }
                return *this;
            }

            /**
            * @brief The value to set the previously selected input.
            * @param[in] value The crisp value to set
            * @return The fisSystem object you invoked << upon.
            */
            fisSystem& operator<<( const real_t& value ) {
                if ( sys.lastTag >= 0 ) {
                    sys.input[ sys.lastTag ].value = value;
                }
                return *this;
            }
    };

    /** @}*/
}

#endif /*QLIBS_FIS*/