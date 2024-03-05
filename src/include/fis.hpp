/*!
 * @file fis.hpp
 * @author J. Camilo Gomez C.
 * @version 1.0.1
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Fuzzy Inference System (FIS) Engine
 **/

#ifndef QLIBS_FIS
#define QLIBS_FIS

#include <include/qlibs_types.hpp>

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {

    /**
    * @brief The Fuzzy Inference System namespace.
    */
    namespace fis {

        /** @addtogroup qfis FIS - Fuzzy Inference System
        * @brief Fuzzy Inference System (FIS) Engine
        * @{
        */

        /**
        * @brief An enum with all the possible values to specify a membership
        * function.
        */
        enum shapeMF {
            custommf = 0,   /*!< Custom user-defined Membership function*/
            trimf,          /*!< Triangular Membership function @c f(a,b,c)*/
            trapmf,         /*!< Trapezoidal Membership function @c f(a,b,c,d)*/
            gbellmf,        /*!< Generalized bell-shaped Membership function @c f(a,b,c)*/
            gaussmf,        /*!< Gaussian  Membership function @c f(s,c)*/
            gauss2mf,       /*!< Gaussian combination Membership function @c f(s1,c1,s2,c2)*/
            sigmf,          /*!< Sigmoidal  Membership function @c f(a,c)*/
            dsigmf,         /*!< Difference between two sigmoidal Membership functions @c f(a1,c1,a2,c2)*/
            psigmf,         /*!< Product of two sigmoidal membership functions @c f(a1,c1,a2,c2)*/
            pimf,           /*!< Pi-shaped membership function @c f(a,b,c,d)*/
            smf,            /*!< S-shaped membership function @c f(a,b)*/
            zmf,            /*!< Z-shaped membership function @c f(a,b)*/
            singletonmf,    /*!< Singleton Membership Function @c f(a)*/
            concavemf,      /*!< Concave Membership Function @c f(i,e)*/
            spikemf,        /*!< Spike Membership Function @c f(w,c)*/
            linsmf,         /*!< Linear s-shaped saturation membership function @c f(a,b)*/
            linzmf,         /*!< Linear z-shaped saturation membership function @c f(a,b)*/
            rectmf,         /*!< Rectangle Membership Function @c f(s,e)*/
            cosmf,          /*!< Cosine Membership Function @c f(c,w)*/
            constantmf,     /*!< Constant membership function @c f(a) [Only for ::Sugeno FIS ]*/
            linearmf,       /*!< Linear membership function @c f(...) [Only for ::Sugeno FIS ]*/
            tlinsmf,        /*!< Tsukamoto s-shaped saturation membership function @c f(a,b) [Only for ::Tsukamoto FIS ]*/
            tlinzmf,        /*!< Tsukamoto linzmf membership function @c f(a,b) [ Only for fis::Tsukamoto FIS ]*/
            tconcavemf,     /*!< Tsukamoto z-shaped saturation membership function @c f(i,e) [Only for ::Tsukamoto FIS ]*/
            tsigmf,         /*!< Tsukamoto Sigmoid membership function @c f(a,c) [ Only for ::Tsukamoto FIS ]*/
            tsmf,           /*!< Tsukamoto S-Shape membership function @c f(a,b) [ Only for ::Tsukamoto FIS ]*/
            tzmf,           /*!< Tsukamoto Z-Shape membership function @c f(a,b) [ Only for ::Tsukamoto FIS ]*/
            trampmf,        /*!< Tsukamoto ramp membership function @c f(a,b) [ Only for ::Tsukamoto FIS ]*/
            /*! @cond  */
            FIS_NUM_MFS        /*!< Number of supported membership functions*/
            /*! @endcond  */
        };

        /**
        * @brief An enum with all the possible de-Fuzzyfication methods.
        */
        enum deFuzzMethod {
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
        enum paramValue {
            FIS_MIN = 0,   /*!< Minimal value*/
            FIS_PROD,      /*!< Product*/
            FIS_MAX,       /*!< Maximum value*/
            FIS_PROBOR,    /*!< Probabilistic OR*/
            FIS_SUM        /*!< Sum*/
        };

        /**
        * @brief An enum with the allowed parameters that can be set on a FIS instance
        */
        enum parameter {
            FIS_Implication,   /*!< Only ::FIS_MIN and ::FIS_PROD supported*/
            FIS_Aggregation,   /*!< Only ::FIS_MAX, ::FIS_PROBOR and ::FIS_SUM supported*/
            FIS_AND,           /*!< Only ::FIS_MIN and ::FIS_PROD supported*/
            FIS_OR,            /*!< Only ::FIS_MAX and ::FIS_PROBOR supported*/
            FIS_EvalPoints     /*!< The number of points for de-fuzzification*/
        };

        /**
        * @brief An enum with the inference system types supported by qlibs::fis
        */
        enum type {
            Mamdani = 0,        /*!< Mamdani inference system. The output of each rule its a fuzzy logic set.*/
            Sugeno,             /*!< Takagi-Sugeno inference system. The output of each rule its a function either linear or constant.*/
            Tsukamoto           /*!< Mamdani inference system. The output of each rule its a fuzzy logic set represented with a monotonic membership function.*/
        };

        /*! @cond  */
        enum deFuzzState {
            FIS_DEFUZZ_INIT,
            FIS_DEFUZZ_COMPUTE,
            FIS_DEFUZZ_END
        };
        /*! @endcond  */

        /*! @cond  */
        class core;
        class instance;
        template<type fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
        class system;
        /*! @endcond  */

        class ioBase {
            protected:
                /*! @cond  */
                real_t min{ -1.0_re };
                real_t max{ 1.0_re };
                real_t value{ 0.0_re };
                /*! @endcond  */
            public:
                virtual ~ioBase() {};
                ioBase() = default;
            friend class instance;
            friend class core;
            template<type fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
            friend class system;
        };

        /**
        * @brief A FIS Input object
        * @details The instance should be initialized using the \ref instance::setupInput() method.
        */
        class input : public ioBase {
            public:
                input() = default;
                virtual ~input() {};
            friend class instance;
            friend class core;
            template<type fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
            friend class system;
        };

        /**
        * @brief A FIS Output object
        * @details The instance should be initialized using the \ref instance::setupOutput() method.
        */
        class output : public ioBase {
            private:
                instance *owner{ nullptr };
                real_t *xag{ nullptr };
                real_t *yag{ nullptr };
                real_t res{ 0.01_re };
                real_t x{ 0.0_re };
                real_t y{ 0.0_re };
                real_t v[ 4 ] = { 0.0_re, 0.0_re, 0.0_re, 0.0_re };
                inline void getNextX( const size_t i ) noexcept
                {
                    x = min + ( ( static_cast<real_t>( i ) + 0.5_re )*res );
                }
            public:
                output() = default;
                virtual ~output() {};
                bool storeAggregatedRegion( real_t *xData,
                                            real_t *yData,
                                            const size_t n ) noexcept;
                template <size_t numberOfPoints>
                bool storeAggregatedRegion( real_t (&xData)[ numberOfPoints ],
                                            real_t (&yData)[ numberOfPoints ] ) noexcept
                {
                    return storeAggregatedRegion( xData, yData, numberOfPoints );
                }
            friend class instance;
            friend class core;
            template<type fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
            friend class system;
        };

        /*! @cond  */
        using mfFunction = real_t (*)( const ioBase * const in,
                                       const real_t *p,
                                       const size_t n );
        /*! @endcond  */

        class core {
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
                virtual ~core() {}
                core() = default;
                static real_t bound( real_t y,
                                     const real_t minValue = 0.0_re,
                                     const real_t maxValue= 1.0_re );
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

                static real_t TriMF( const ioBase * const in,
                                     const real_t *p,
                                     const size_t n );
                static real_t TrapMF( const ioBase * const in,
                                      const real_t *p,
                                      const size_t n );
                static real_t GBellMF( const ioBase * const in,
                                       const real_t *p,
                                       const size_t n );
                static real_t GaussMF( const ioBase * const in,
                                       const real_t *p,
                                       const size_t n );
                static real_t Gauss2MF( const ioBase * const in,
                                        const real_t *p,
                                        const size_t n );
                static real_t SigMF( const ioBase * const in,
                                     const real_t *p,
                                     const size_t n );
                static real_t TSigMF( const ioBase * const in,
                                      const real_t *p,
                                      const size_t n );
                static real_t DSigMF( const ioBase * const in,
                                      const real_t *p,
                                      const size_t n );
                static real_t PSigMF( const ioBase * const in,
                                      const real_t *p,
                                      const size_t n );
                static real_t SMF( const ioBase * const in,
                                   const real_t *p,
                                   const size_t n );
                static real_t TSMF( const ioBase * const in,
                                    const real_t *p,
                                    const size_t n );
                static real_t ZMF( const ioBase * const in,
                                   const real_t *p,
                                   const size_t n );
                static real_t LinSMF( const ioBase * const in,
                                      const real_t *p,
                                      const size_t n );
                static real_t LinZMF( const ioBase * const in,
                                      const real_t *p,
                                      const size_t n );
                static real_t TZMF( const ioBase * const in,
                                    const real_t *p,
                                    const size_t n );
                static real_t PiMF( const ioBase * const in,
                                    const real_t *p,
                                    const size_t n );
                static real_t SingletonMF( const ioBase * const in,
                                           const real_t *p,
                                           const size_t n );
                static real_t ConcaveMF( const ioBase * const in,
                                         const real_t *p,
                                         const size_t n );
                static real_t TConcaveMF( const ioBase * const in,
                                          const real_t *p,
                                          const size_t n );
                static real_t SpikeMF( const ioBase * const in,
                                       const real_t *p,
                                       const size_t n );
                static real_t TLinSMF( const ioBase * const in,
                                       const real_t *p,
                                       const size_t n );
                static real_t TLinZMF( const ioBase * const in,
                                       const real_t *p,
                                       const size_t n );
                static real_t TRampMF( const ioBase * const in,
                                       const real_t *p,
                                       const size_t n );
                static real_t RectangleMF( const ioBase * const in,
                                           const real_t *p,
                                           const size_t n );
                static real_t CosineMF( const ioBase * const in,
                                        const real_t *p,
                                        const size_t n );
                static real_t ConstantMF( const ioBase * const in,
                                          const real_t *p,
                                          const size_t n );
                static real_t LinearMF( const ioBase * const in,
                                        const real_t *p,
                                        const size_t n );

                static real_t deFuzzCentroid( output * const o,
                                            const deFuzzState stage );
                static real_t deFuzzBisector( output * const o,
                                            const deFuzzState stage );
                static real_t deFuzzLOM( output * const o,
                                        const deFuzzState stage );
                static real_t deFuzzSOM( output * const o,
                                        const deFuzzState stage );
                static real_t deFuzzMOM( output * const o,
                                        const deFuzzState stage );
                static real_t deFuzzWtAverage( output * const o,
                                            const deFuzzState stage );
                static real_t deFuzzWtSum( output * const o,
                                        const deFuzzState stage );
                /*! @endcond  */
        };

        /**
        * @brief A FIS Membership Function
        * @details The instance should be initialized using the \ref instance::setupInputMF()
        * or \ref instance::setupOutputMF() methods.
        */
        class mf {
            private:
                mfFunction shape{ nullptr };
                const real_t *points{ nullptr };
                real_t fx{ 0.0_re };
                real_t h{ 0.0_re };
                size_t index{ 0U };
                inline size_t getIndex( void ) const
                {
                    return index;
                }
                inline real_t membership( const ioBase * const io,
                                          size_t n = 1U )
                {
                    fx = ( nullptr != shape ) ? h*shape( io, points, n ) : 0.0_re;
                    return fx;
                }
            public:
                virtual ~mf() {};
                mf() = default;
                friend class instance;
        };

        /*! @cond  */
        #define FIS_RULE_ITEM_SIZE  1
        /*! @endcond  */

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
            * static const fis::rules Rules[] = {
            *     FIS_RULES_BEGIN
            *         IF service IS service_poor OR food IS food_rancid THEN tip IS tip_cheap END
            *         IF service IS service_good THEN tip IS tip_average END
            *         IF service IS service_excellent OR food IS food_delicious THEN tip IS tip_generous END
            *     FIS_RULES_END
            * };
            * @endcode
            */
            using rules = int8_t;
            /*! @cond  */
            constexpr fis::rules FIS_RULES_MIN_VALUE = INT8_MIN;
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
            * static const fis::rules Rules[] = {
            *     FIS_RULES_BEGIN
            *         IF service IS service_poor OR food IS food_rancid THEN tip IS tip_cheap END
            *         IF service IS service_good THEN tip IS tip_average END
            *         IF service IS service_excellent OR food IS food_delicious THEN tip IS tip_generous END
            *     FIS_RULES_END
            * };
            * @endcode
            */
            using rules = int16_t;
            /*! @cond  */
            constexpr fis::rules FIS_RULES_MIN_VALUE = INT16_MIN;
            /*! @endcond  */
        #endif

        /**
        * @brief Used to define an enum of fis tags
        * @details Tags for I/O and set/membership functions should be defined
        * within enums of this type
        * Example:
        * @code{.c}
        * enum : fis::tag { tag1, tag2, tag3 };
        * @endcode
        */
        using tag = rules;
        /*cstat -MISRAC++2008-0-1-4_b*/
        /*! @cond  */
        constexpr fis::rules Q_FIS_RULES_BEGIN    = FIS_RULES_MIN_VALUE;
        constexpr fis::rules Q_FIS_RULES_END      = FIS_RULES_MIN_VALUE + 1;
        constexpr fis::rules Q_FIS_AND            = FIS_RULES_MIN_VALUE + 2;
        constexpr fis::rules Q_FIS_OR             = FIS_RULES_MIN_VALUE + 3;
        constexpr fis::rules Q_FIS_THEN           = FIS_RULES_MIN_VALUE + 4;
        /*cstat +MISRAC++2008-0-1-4_b*/
        #define Q_FIS_IF_STATEMENT      ,
        #define Q_FIS_AND_STATEMENT     +1),fis::Q_FIS_AND,
        #define Q_FIS_OR_STATEMENT      +1),fis::Q_FIS_OR,
        #define Q_FIS_THEN_STATEMENT    +1),fis::Q_FIS_THEN,
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
        * static const fis::rules Rules[] = {
        *     FIS_RULES_BEGIN
        *
        *     FIS_RULES_END
        * };
        * @endcode
        */
        #define FIS_RULES_BEGIN         fis::Q_FIS_RULES_BEGIN
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
        * static const fis::rules Rules[] = {
        *     FIS_RULES_BEGIN
        *
        *     FIS_RULES_END
        * };
        * @endcode
        */
        #define FIS_RULES_END           ,fis::Q_FIS_RULES_END
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


    using deFuzzFunction = real_t (*)( output * const o, const deFuzzState stage );
    using fuzzyOperator = real_t (*)( const real_t a, const real_t b );

    /**
    * @brief A FIS(Fuzzy Inference System) object
    * @details The instance should be initialized using the \ref instance::setup() method.
    */
     class instance: public core, private nonCopyable {
        using methods_fcn = real_t (*)( const real_t a, const real_t b );

        template<type fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
        friend class system;

        private:
            input *xInput{ nullptr };
            output *xOutput{ nullptr };
            mf *inMF{ nullptr };
            mf *outMF{ nullptr };
            methods_fcn andOp{ &Min };
            methods_fcn orOp{ &Max };
            methods_fcn implicate{ &Min };
            methods_fcn aggregate{ &Max };
            size_t (instance::*inferenceState)( size_t i );
            size_t (instance::*aggregationState)( size_t i );
            deFuzzFunction deFuzz{ &deFuzzCentroid };
            real_t *ruleWeight{ nullptr };
            real_t *wi{ nullptr };
            const rules *xRules{ nullptr };
            size_t rules_cols{ 0U};
            size_t nInputs{ 0 };
            size_t nOutputs{ 0 };
            size_t nMFInputs{ 0U };
            size_t nMFOutputs{ 0U };
            size_t nPoints{ 100 };
            size_t nRules{ 0U };
            size_t ruleCount{ 0U };
            real_t rStrength{ 0.0_re };
            rules lastConnector;
            type xType{ Mamdani };
            bool setMF( mf *m,
                        const tag io,
                        const tag mf,
                        const shapeMF s,
                        mfFunction customMf,
                        const real_t *cp,
                        const real_t h ) noexcept;
            static real_t parseFuzzValue( mf * const mfIO,
                                          rules index ) noexcept;
            fuzzyOperator getFuzzOperator( void ) noexcept;
            size_t inferenceAntecedent( size_t i ) noexcept;
            size_t inferenceReachEnd( size_t i ) noexcept;
            size_t aggregationFindConsequent( size_t i ) noexcept;
            size_t inferenceConsequent( size_t i ) noexcept;
            void fuzzyAggregate( void ) noexcept;
            static const size_t INFERENCE_ERROR;
            tag lastTag{ -1 };
        public:
            instance() = default;
            virtual ~instance() {}

            /**
            * @brief Setup and initialize the FIS instance.
            * @note Default configuration : AND = Min, OR = Max, Implication = Min
            * Aggregation = Max, EvalPoints = 100
            * @param[in] t Type of inference ::Mamdani, ::Sugeno or ::Tsukamoto.
            * @param[in] inputs An array with all the system inputs as fis::input
            * objects.
            * @param[in] ni The number of elements in the @a inputs array.
            * @param[in] outputs An array with all the system outputs as fis::output
            * objects.
            * @param[in] no The number of elements in the @a outputs array.
            * @param[in] mf_inputs An array with all the membership functions related to
            * the inputs. This should be an array of fis::mf objects.
            * @param[in] nmi The number of elements in the @a mf_inputs array.
            * @param[in] mf_outputs An array with all the membership functions related to
            * the outputs. This should be an array of fis::mf objects.
            * @param[in] nmo The number of elements in the @a mf_outputs array.
            * @param[in] r The rules set.
            * @param[in] n Number of rules
            * @param[in] rWeights An array of size @a n were the rule strengths will be stored.
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( const type t,
                 input * const inputs,
                 const size_t ni,
                 output * const outputs,
                 const size_t no,
                 mf * const mf_inputs,
                 const size_t nmi,
                 mf * const mf_outputs,
                 const size_t nmo,
                 const rules * const r,
                 const size_t n,
                 real_t *rWeights = nullptr ) noexcept;

            /**
            * @brief Setup and initialize the FIS instance.
            * @note Default configuration : AND = Min, OR = Max, Implication = Min
            * Aggregation = Max, EvalPoints = 100
            * @param[in] t Type of inference ::Mamdani, ::Sugeno or ::Tsukamoto.
            * @param[in] inputs An array with all the system inputs as fis::input
            * objects.
            * @param[in] outputs An array with all the system outputs as fis::output
            * objects.
            * @param[in] mf_inputs An array with all the membership functions related to
            * the inputs. This should be an array of fis::mf objects.
            * @param[in] mf_outputs An array with all the membership functions related to
            * the outputs. This should be an array of fis::mf objects.
            * @param[in] r The rules set.
            * @param[in] rWeights An array of size @a n were the rule strengths will be stored.
            * @return @c true on success, otherwise return @c false.
            */
            template <size_t numberInputs, size_t numberOutputs, size_t numberMFinputs, size_t numberMFOutputs, size_t numberRules>
            bool setup( const type t,
                 input (&inputs)[ numberInputs ],
                 output (&outputs)[ numberOutputs ],
                 mf (&mf_inputs)[ numberMFinputs ],
                 mf (&mf_outputs)[ numberMFOutputs ],
                 const rules * const r,
                 real_t (&rWeights)[ numberRules ] ) noexcept
            {
                return setup( t, inputs, numberInputs, outputs, numberOutputs, mf_inputs, numberMFinputs, mf_outputs, numberMFOutputs, r, numberRules, rWeights );
            }

            /**
            * @brief Setup the input with the specified tag and set limits for it
            * @param[in] t The input tag
            * @param[in] Min Minimum allowed value for this input
            * @param[in] Max Max allowed value for this input
            * @return @c true on success, otherwise return @c false.
            */
            bool setupInput( const tag t,
                             const real_t Min,
                             const real_t Max ) noexcept;
            /**
            * @brief Setup the output with the specified tag and set limits for it
            * @param[in] t The output tag
            * @param[in] Min Minimum allowed value for this output
            * @param[in] Max Max allowed value for this output
            * @return @c true on success, otherwise return @c false.
            */
            bool setupOutput( const tag t,
                              const real_t Min,
                              const real_t Max ) noexcept;

            /**
            * @brief Setup the input tag and points for the specified membership
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
            bool setupInputMF( const tag io,
                               const tag mf,
                               const shapeMF s,
                               const real_t *cp,
                               const real_t h = 1.0_re ) noexcept
            {
                return ( ( nullptr != inMF ) && ( nMFInputs > 0U ) ) ? setMF( inMF, io, mf, s, nullptr, cp, h ) : false;
            }

            /**
            * @brief Setup the input tag and points for the specified membership
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
            bool setupInputMF( const tag io,
                               const tag mf,
                               mfFunction customMfs,
                               const real_t *cp,
                               const real_t h = 1.0_re ) noexcept
            {
                return ( ( nullptr != inMF ) && ( nMFInputs > 0U ) ) ? setMF( inMF, io, mf, custommf, customMfs, cp, h ) : false;
            }

            /**
            * @brief Setup the output tag and points for the specified membership
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
            bool setupOutputMF( const tag io,
                                const tag mf,
                                const shapeMF s,
                                const real_t *cp,
                                const real_t h = 1.0_re ) noexcept
            {
                return ( ( nullptr != outMF ) && ( nMFOutputs > 0U ) ) ? setMF( outMF, io, mf, s, nullptr, cp, ( ( Mamdani == xType ) ? h : 1.0_re ) ) : false;
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
            bool setupOutputMF( const tag io,
                                const tag mf,
                                mfFunction customMfs,
                                const real_t *cp,
                                const real_t h = 1.0_re ) noexcept
            {
                return ( ( nullptr != outMF ) && ( nMFOutputs > 0U ) ) ? setMF( outMF, io, mf, custommf, customMfs, cp, ( ( Mamdani == xType ) ? h : 1.0_re ) ) : false;
            }

            /**
            * @brief Set a crisp value of the input with the specified tag.
            * @param[in] t The input tag
            * @param[in] value The crisp value to set
            * @return @c true on success, otherwise return @c false.
            */
            bool setInput( const tag t,
                           const real_t value ) noexcept;

            /**
            * @brief Get the de-fuzzified crisp value from the output with the
            * specified tag.
            * @param[in] t The output tag
            * @param[out] value The variable where the output will be stored
            * @return The requested de-fuzzified crisp value.
            */
            bool getOutput( const tag t,
                            real_t &value ) const noexcept;

            /**
            * @brief Set parameters of the FIS instance.
            * @param[in] p The requested parameter to change/set.
            * @param[in] x The value of the parameter to set.
            * @return @c true on success, otherwise return @c false.
            */
            bool setParameter( const parameter p,
                               const paramValue x ) noexcept;

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
            bool setDeFuzzMethod( deFuzzMethod m ) noexcept;

            /**
            * @brief Perform the fuzzification operation over the crisp inputs on the
            * requested FIS object
            * @pre I/Os and fuzzy sets must be previously initialized by instance::setupInput(),
            * instance::setupOutput(), instance::setupInputMF(), instance::setupOutputMF and instance::setup() respectively.
            * @return @c true on success, otherwise return @c false.
            */
            bool fuzzify( void ) noexcept;

            /**
            * @brief Perform the de-Fuzzification operation to compute the crisp outputs.
            * @pre The instance should have already invoked the inference process
            * successfully with fis::inference()
            * @note By default, this method, uses the Centroid method on
            * ::Mamdani type FIS and weight-average on ::Sugeno type FIS. To change
            * the default settings use the instance::setDeFuzzMethod() function.
            * @return @c true on success, otherwise return @c false.
            */
            bool deFuzzify( void ) noexcept;

            /**
            * @brief Perform the inference process on the FIS object
            * @pre The instance should have already invoked the fuzzification operation
            * successfully with instance::fuzzify()
            * @return @c true on success, otherwise return @c false.
            */
            bool inference( void ) noexcept;

            /**
            * @brief Set weights to the rules of the inference system.
            * @pre I/Os and fuzzy sets must be previously initialized by instance::setupInput(),
            * instance::setupOutput(), instance::setInputMF(), instance::setOutputMF and instance::setup() respectively.
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
            real_t operator[]( tag outTag ) const
            {
                /*cstat -CERT-STR34-C*/
                return ( static_cast<size_t>( outTag ) < nOutputs ) ? xOutput[ outTag ].value :  xInput[ nOutputs -1 ].value;
                /*cstat +CERT-STR34-C*/
            }

            /**
            * @brief Select the input to set using with the specified tag.
            * @param[in] t The input tag
            * @return The fis object you invoked << upon.
            */
            instance& operator<<( const tag& t ) {
                if ( static_cast<size_t>( t ) < nInputs ) {
                    lastTag = t;
                }
                return *this;
            }

            /**
            * @brief The value to set the previously selected input.
            * @param[in] value The crisp value to set
            * @return The fis object you invoked << upon.
            */
            instance& operator<<( const int& value ) {
                if ( lastTag >= 0 ) {
                    xInput[ lastTag ].value = static_cast<real_t>(value);
                }
                return *this;
            }

            /**
            * @brief The value to set the previously selected input.
            * @param[in] value The crisp value to set
            * @return The fis object you invoked << upon.
            */
            instance& operator<<( const real_t& value ) {
                if ( lastTag >= 0 ) {
                    xInput[ lastTag ].value = value;
                }
                return *this;
            }

        friend class core;
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
    template<type fType, size_t numberOfInputs, size_t numberOfOutputs, size_t numberOfInputSets, size_t numberOfOutputSets, size_t numberOfRules>
    class system {
        private:
            instance sys;
            const rules *xRules;
            input inputs[ numberOfInputs ];
            output outputs[ numberOfOutputs ];
            mf MFin[ numberOfInputSets ], MFout[ numberOfOutputSets ];
            real_t ruleStrength[ numberOfRules ];
        public:
            constexpr system( const rules *r ) : xRules( r ) {};

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
                                  xRules, numberOfRules,
                                  ruleStrength );
            }

            /**
            * @brief Setup the input with the specified tag and set limits for it
            * @param[in] t The input tag
            * @param[in] Min Minimum allowed value for this input
            * @param[in] Max Max allowed value for this input
            * @return @c true on success, otherwise return @c false.
            */
            inline bool setupInput( const tag t,
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
            inline bool setupOutput( const tag t,
                                     const real_t Min,
                                     const real_t Max ) noexcept
            {
                return sys.setupOutput( t, Min, Max );
            }

            /**
            * @brief Setup the input tag and points for the specified membership
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
            inline bool setupInputMF( const tag io,
                                      const tag mf,
                                      const shapeMF s,
                                      const real_t *cp,
                                      const real_t h = 1.0_re ) noexcept
            {
                return sys.setupInputMF( io, mf, s, cp, h );
            }

            /**
            * @brief Setup the input tag and points for the specified membership
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
            inline bool setupInputMF( const tag io,
                                      const tag mf,
                                      mfFunction customMfs,
                                      const real_t *cp,
                                      const real_t h = 1.0_re ) noexcept
            {
                return sys.setupInputMF( io, mf, customMfs, cp, h );
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
            inline bool setupOutputMF( const tag io,
                                       const tag mf,
                                       const shapeMF s,
                                       const real_t *cp,
                                       const real_t h = 1.0_re ) noexcept
            {
                return sys.setupOutputMF( io, mf, s, cp, h );
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
            inline bool setupOutputMF( const tag io,
                                       const tag mf,
                                       mfFunction customMfs,
                                       const real_t *cp,
                                       const real_t h = 1.0_re ) noexcept
            {
                return sys.setupOutputMF( io, mf, customMfs, cp, h );
            }

            /**
            * @brief Perform the fuzzification operation over the crisp inputs on the
            * requested FIS object
            * @pre I/Os and fuzzy sets must be previously initialized by system::setupInput(),
            * system::setupOutput(), system::setupInputMF(), system::setupOutputMF
            * and system::setup() respectively.
            * @return @c true on success, otherwise return @c false.
            */
            inline bool fuzzify( void ) noexcept
            {
                return sys.fuzzify();
            }

            /**
            * @brief Perform the de-Fuzzification operation to compute the crisp outputs.
            * @pre The instance should have already invoked the inference process
            * successfully with system::inference()
            * @note By default, this method, uses the Centroid method on
            * ::Mamdani type FIS and weight-average on ::Sugeno type FIS. To change
            * the default settings use the system::setDeFuzzMethod() function.
            * @return @c true on success, otherwise return @c false.
            */
            inline bool deFuzzify( void ) noexcept
            {
                return sys.deFuzzify();
            }

            /**
            * @brief Perform the inference process on the FIS object
            * @pre The instance should have already invoked the fuzzification operation
            * successfully with system::fuzzify()
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
            inline bool setInput( const tag t,
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
            inline bool getOutput( const tag t,
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
            inline bool setParameter( const parameter p,
                                      const paramValue x ) noexcept
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
            inline bool setDeFuzzMethod( deFuzzMethod m ) noexcept
            {
                return sys.setDeFuzzMethod( m );
            }

            /**
            * @brief Set weights to the rules of the inference system.
            * @pre I/Os and fuzzy sets must be previously initialized by
            * system::setupInput(), system::setupOutput(),
            * system::setInputMF(), system::setOutputMF and
            * system::setup() respectively.
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
            real_t operator[]( tag outTag ) const
            {
                /*cstat -CERT-STR34-C*/
                return ( static_cast<size_t>( outTag ) < sys.nOutputs ) ? sys.xOutput[ outTag ].value :  sys.xInput[ sys.nOutputs -1 ].value;
                /*cstat +CERT-STR34-C*/
            }

            /**
            * @brief Select the input to set using with the specified tag.
            * @param[in] t The input tag
            * @return The fis::system object you invoked << upon.
            */
            system& operator<<( const tag& t ) {
                if ( static_cast<size_t>( t ) < sys.nInputs ) {
                    sys.lastTag = t;
                }
                return *this;
            }

            /**
            * @brief The value to set the previously selected input.
            * @param[in] value The crisp value to set
            * @return The fis::system object you invoked << upon.
            */
            system& operator<<( const int& value ) {
                if ( sys.lastTag >= 0 ) {
                    sys.xInput[ sys.lastTag ].value = static_cast<real_t>( value );
                }
                return *this;
            }

            /**
            * @brief The value to set the previously selected input.
            * @param[in] value The crisp value to set
            * @return The fis::system object you invoked << upon.
            */
            system& operator<<( const real_t& value ) {
                if ( sys.lastTag >= 0 ) {
                    sys.xInput[ sys.lastTag ].value = value;
                }
                return *this;
            }
    };
        /** @}*/
    }

}

#endif /*QLIBS_FIS*/