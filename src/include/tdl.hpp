/*!
 * @file tdl.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Helper class that implements a Tapped Delay Line (TDL). A TDL is a
 * delay line that provides access to its contents at arbitrary intermediate
 * delay length values.
 * This class runs in constant time O(1), so it becomes useful when you need to
 * work with long delayed lines.
 **/

#ifndef QLIBS_TDL
#define QLIBS_TDL

#include "include/types.hpp" 

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {
    /** @addtogroup  qtdl Tapped Delay Line
    * @brief An implementation of the Tapped Delay Line (TDL) structure in O(1) 
    *  @{
    */


    /**
    * @brief A Tapped Delay Line (TDL) object
    * @details The instance should be initialized using the td::setup() method.
    */
    class tdl : private nonCopyable {
        protected:
            /*! @cond  */
            real_t *head{ nullptr };
            real_t *tail{ nullptr };
            real_t *rd{ nullptr };
            real_t *wr{ nullptr };
            size_t itemCount{ 0U };
            const real_t undefined{ nan("" ) }; // skipcq: CXX-W2010

            void insertNewest( const real_t sample ) noexcept;
            void removeOldest( void ) noexcept;
            /*! @endcond  */
        public:
            virtual ~tdl() {}
            tdl() = default;

            /**
            * @brief Constructor for the Tapped Delay Line (TDL) instance.
            * @param[in] area An array of size @a n where delays will be stored
            * @param[in] n The number of elements on @a area.
            * @param[in] initVal The value with which all TDL delays will be initialized
            */
            tdl( real_t * const area,
                 const size_t n,
                 const real_t initVal = 0.0 )
            {
                setup( area, n, initVal );
            }

            /**
            * @brief Constructor for the Tapped Delay Line (TDL) instance.
            * @param[in] area An array where delays will be stored
            * @param[in] initVal The value with which all TDL delays will be initialized
            */
            template <size_t numberOfDelays>
            tdl( real_t (&area)[ numberOfDelays ],
                 const real_t initVal = 0.0 ) noexcept
            {
                setup( area, numberOfDelays, initVal );
            }

            /**
            * @brief Setup and initialize a Tapped Delay Line (TDL) instance by setting
            * the default optimal parameters.
            * @param[in] area An array of size @a n where delays will be stored
            * @param[in] n The number of elements on @a area.
            * @param[in] initVal The value with which all TDL delays will be initialized
            * @return none
            */
            void setup( real_t * const area,
                        const size_t n,
                        const real_t initVal = 0.0 ) noexcept;

            /**
            * @brief Setup and initialize a Tapped Delay Line (TDL) instance by setting
            * the default optimal parameters.
            * @param[in] area The array where delays will be stored
            * @param[in] initVal The value with which all TDL delays will be initialized
            * @return none
            */
            template <size_t numberOfDelays>
            void setup( real_t (&area)[ numberOfDelays ],
                        const real_t initVal = 0.0 ) noexcept
            {
                setup( area, numberOfDelays, initVal );
            }

            /**
            * @brief Clears all delays from the TDL and sets them to the specified value
            * @param[in] initVal The value with which all TDL delays will be initialized
            * @return none
            */
            void flush( const real_t initVal = 0.0 ) noexcept;

            /**
            * @brief Get the oldest sample from the TDL x(k-n)
            * @return The requested value from the TDL
            */
            real_t getOldest( void ) const noexcept;

            /**
            * @brief Get the most recent sample from the TDL x(k)
            * @return The requested value from the TDL
            */
            real_t getRecent( void ) const noexcept;

            /**
            * @brief Get the specified delayed sample from the TDL x(k-i)
            * @param[in] i The requested delay index
            * @return The requested value from the TDL
            */
            real_t getAtIndex( const size_t i ) const noexcept;

            /**
            * @brief Insert a new sample to the TDL removing the oldest sample
            * @param[in] sample The new sample to insert.
            * @return none
            */
            void insertSample( const real_t sample ) noexcept;

            /**
            * @brief Get the specified delayed sample from the TDL x(k-i)
            * @param[in] index The requested delay index
            * @return The requested value from the TDL
            */
            const real_t& operator[]( int index ) noexcept;

            /**
            * @brief Insert a new sample to the TDL removing the oldest sample
            * @param[in] sample The new sample to insert.
            * @return none
            */
            void operator()( const real_t sample ) noexcept
            {
                insertSample( sample );
            }

            /**
            * @brief Check if the TDL has been initialized.
            * @return @c true if instance has been initialized
            */
            bool isInitialized( void ) const {
                return ( nullptr != head );
            }
    };

    /** @}*/
}


#endif /*QLIBS_TDL*/