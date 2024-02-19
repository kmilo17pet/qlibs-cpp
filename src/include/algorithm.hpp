#ifndef QLIBS_ALGORITHM
#define QLIBS_ALGORITHM

#include <include/qlibs_types.hpp>

/*!
 * @file algorithm.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief A basic implementation of basic algorithms for raw-arrays
 **/


/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {

    /**
    * @brief The generic namespace.
    */
    namespace algorithm {
        /** @addtogroup  qalgorithm Basic-algorithms for raw-arrays
        * @brief A basic implementation of basic algorithms for raw-arrays
        * without recursion and dynamic memory allocation
        * @{
        */

        /**
        * @brief Exchanges the values of a and b.
        * @param[in,out] x Object to be swapped.
        * @param[in,out] y Object to be swapped.
        * @return none.
        */
        template <typename T>
        void swap( T& x, T& y ) noexcept
        {
            T tmp = x;
            x = y;
            y = tmp;
        }

         /** @cond */
        namespace impl {
            template<typename T, size_t N>
            class sort_stack {
                private:
                    T dat[ N ];
                    size_t topIndex;
                public:
                    sort_stack() : topIndex( 0 ) {}
                    bool empty( void ) const
                    {
                        return 0 == topIndex;
                    }
                    void push( const T& value )
                    {
                        if ( topIndex < N ) {
                            dat[ topIndex++ ] = value;
                        }
                    }
                    void pop( void )
                    {
                        if ( topIndex > 0 ) {
                            --topIndex;
                        }
                    }
                    T& top( void )
                    {
                        return dat[ topIndex - 1 ];
                    }
            };
        }
         /** @endcond */

        /**
        * @brief Sorts the given array in the range [first,last) into ascending
        * order.
        * @note The elements are compared using operator<
        * @remark This algorithm uses a non-recursive variant of the quicksort
        * algorithm.
        * @param[in,out] array The array to be sorted.
        * @param[in] first Initial position of the portion to be sorted
        * @param[in] last Final position of the portion to be sorted
        * @return none.
        */
        template<typename T, size_t n>
        void sort( T ( &array )[ n ],
                       size_t first = 0U,
                       size_t last = n - 1U )
        {
            struct pair {
                int first;
                int second;
            };
            algorithm::impl::sort_stack<pair, n> stack;
            int start = static_cast<int>( first );
            int end = static_cast<int>( last );

            stack.push( { start, end } );
            while ( !stack.empty() ) {
                pair indices = stack.top();
                stack.pop();
                start = indices.first;
                end = indices.second;
                int pivotIndex = start;
                T pivotValue = array[ end ];

                for ( int i = start; i < end; ++i ) {
                    if ( array[ i ] < pivotValue ) {
                        algorithm::swap( array[ i ], array[ pivotIndex ] );
                        ++pivotIndex;
                    }
                }
                algorithm::swap( array[ pivotIndex ],  array[ end ] );
                if ( pivotIndex - 1 > start ) {
                    stack.push( { start, pivotIndex - 1 } );
                }
                if ( pivotIndex + 1 < end ) {
                    stack.push( { pivotIndex + 1, end } );
                }
            }
        }

        /**
        * @brief Reverse the given array. Operation takes place on the portion
        * of the array that starts at position @a init to position @a end.
        * @param[in,out] array The array to reverse.
        * @param[in] init Position of the first element.
        * @param[in] end Position of the last element.
        * @return none.
        */

        template<typename T, size_t n>
        void reverse( T ( &array )[ n ],
                      const size_t init = 0U,
                      const size_t end = n - 1U ) noexcept
        {
            if ( end > init ) {
                size_t s = init, e = end;

                while( s < e ) {
                    algorithm::swap( array[ s ], array[ e ] );
                    ++s;
                    --e;
                }
            }
        }

        /**
        * @brief Rotates @a k elements of the array pointed. Rotation direction
        * is determined by the sign of @a k, the means a positive value performs
        *  a right-rotation and a negative value a left-rotation.
        * @param[in,out] array The array to rotate.
        * @param[in] k Positions to rotate.
        * @return none.
        */

        template<typename T, size_t n>
        void rotate( T ( &array )[ n ],
                         const int k = 1 ) noexcept
        {
            if ( 0 != k ) {
                size_t r;
                if ( k > 0 ) {
                    r = static_cast<size_t>( k );
                    r %= n;
                    algorithm::reverse( array, n - r, n - 1U );
                    algorithm::reverse( array, 0U, n - r - 1U );
                    algorithm::reverse( array, 0U, n - 1U );
                }
                else {
                    /*cstat -MISRAC++2008-5-0-9*/
                    r = static_cast<size_t>( -k );
                    /*cstat +MISRAC++2008-5-0-9*/
                    r %= n;
                    algorithm::reverse( array, 0U, r - 1U );
                    algorithm::reverse( array, r, n - 1U );
                    algorithm::reverse( array, 0U, n - 1U );
                }
            }
        }

        /**
         * @brief Sets all elements of an array in the range  [first,last) to a
         * specific value.
         * @param[in,out] array The array to set.
         * @param[in] value The value to set all elements to.
         * @param[in] first Initial position of the portion to search
         * @param[in] last Final position of the portion to search
         */
        template<typename T, size_t n>
        inline void set( T ( &array )[ n ],
                         const T value,
                         const size_t first = 0U,
                         const size_t last = n - 1U ) noexcept
        {
            for ( size_t i = first ; i < last; ++i ) {
                array[ i ] = value;
            }
        }

        /**
        * @brief Performs a linear search over the raw-array in the
        * range [first,last) that matches the @a key.
        * @note The elements are compared using operator '=='
        * @param[in] key The object that serves as key for
        * the search.
        * @param[in] array The array where the search is performed
        * @param[in] first Initial position of the portion to search
        * @param[in] last Final position of the portion to search
        * @return This function returns a pointer to an entry in the array that
        * matches the search key. If key is not found, a @c nullptr pointer is
        * returned.
        */
        template<typename T, size_t n>
        inline T* lSearch( const T key,
                           T ( &array )[ n ],
                           const size_t first = 0U,
                           const size_t last = n - 1U
                           ) noexcept
        {
            T* found = nullptr;

            for ( size_t i = first; i < last; ++i ) {
                if ( array[ i ] == key ) {
                    found = &array[ i ];
                    break;
                }
            }
            return found;
        }

        /**
        * @brief Performs a binary search over the raw-array in the
        * range [first,last) that matches the @a key..
        * The array contents should be sorted in ascending order.
        * @note The elements are compared using operator '<' and '=='
        * @param[in] key The object that serves as key for
        * the search.
        * @param[in] array The array where the search is performed
        * @param[in] first Initial position of the portion to search
        * @param[in] last Final position of the portion to search
        * @return This function returns a pointer to an entry in the array that
        * matches the search key. If key is not found, a @c nullptr pointer is returned.
        */
        template<typename T, size_t n>
        inline T* bSearch( const T key,
                           T ( &array )[ n ],
                           const size_t first = 0U,
                           const size_t last = n - 1U
                           ) noexcept
        {
            T* found = nullptr;
            int left = static_cast<int>( first );
            int right = static_cast<int>( last );

            while ( left <= right ) {
                int mid = left + ( right - left )/2;

                if ( array[ mid ] == key ) {
                    found =  &array[ mid ];
                    break;
                }
                if ( array[ mid ] < key ) {
                    left = mid + 1;
                }
                else {
                    right = mid - 1;
                }
            }
            return found;
        }

        /** @}*/
    }
}


#endif /*QLIBS_ALGORITHM*/