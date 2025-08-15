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

#include <include/qlibs_types.hpp>

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

            void insertNewest( const real_t sample ) noexcept;
            void removeOldest( void ) noexcept;
            /*! @endcond  */
        public:
            virtual ~tdl() noexcept = default;
            tdl() = default;

            /**
            * @brief Constructor for the Tapped Delay Line (TDL) instance.
            * @param[in] area An array of size @a n where delays will be stored
            * @param[in] n The number of elements on @a area.
            * @param[in] initVal The value with which all TDL delays will be initialized
            */
            tdl( real_t * const area,
                 const size_t n,
                 const real_t initVal = 0.0_re ) noexcept
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
                 const real_t initVal = 0.0_re ) noexcept
            {
                setup( area, numberOfDelays, initVal );
            }

            /**
            * @brief Setup and initialize a Tapped Delay Line (TDL) instance by setting
            * the default optimal parameters.
            * @param[in] area An array of size @a n where delays will be stored
            * @param[in] n The number of elements on @a area.
            * @param[in] initVal The value with which all TDL delays will be initialized
            */
            void setup( real_t * const area,
                        const size_t n,
                        const real_t initVal = 0.0_re ) noexcept;

            /**
            * @brief Setup and initialize a Tapped Delay Line (TDL) instance by setting
            * the default optimal parameters.
            * @param[in] area The array where delays will be stored
            * @param[in] initVal The value with which all TDL delays will be initialized
            */
            template <size_t numberOfDelays>
            void setup( real_t (&area)[ numberOfDelays ],
                        const real_t initVal = 0.0_re ) noexcept
            {
                setup( area, numberOfDelays, initVal );
            }

            /**
            * @brief Clears all delays from the TDL and sets them to the specified value
            * @param[in] initVal The value with which all TDL delays will be initialized
            */
            void flush( const real_t initVal = 0.0_re ) noexcept;

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
            */
            void insertSample( const real_t sample ) noexcept;

            /**
            * @brief Get the specified delayed sample from the TDL x(k-i)
            * @param[in] index The requested delay index
            * @return The requested value from the TDL
            */
            real_t operator[]( int index ) noexcept;

            /**
            * @brief Insert a new sample to the TDL removing the oldest sample
            * @param[in] sample The new sample to insert.
            */
            void operator()( const real_t sample ) noexcept
            {
                insertSample( sample );
            }

            /**
            * @brief Check if the TDL has been initialized.
            * @return @c true if instance has been initialized
            */
            bool isInitialized( void ) const noexcept {
                return ( nullptr != head );
            }

            /**
            * @brief Check if the TDL has been initialized.
            * @return @c true if instance has been initialized
            */
            explicit operator bool() const noexcept {
                return isInitialized(); // controls truthiness
            }
    };

    /** @}*/
}


#endif /*QLIBS_TDL*/