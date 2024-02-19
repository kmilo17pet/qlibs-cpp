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
        * @brief Reverses the order of the elements in the range [first,last).
        * @param[in,out] array The array to reverse.
         * @param[in] first Initial position of the portion to reverse
         * @param[in] last Final position of the portion to reverse
        * @return none.
        */

        template<typename T, size_t n>
        void reverse( T ( &array )[ n ],
                      const size_t first = 0U,
                      const size_t last = n - 1U ) noexcept
        {
            if ( last > first ) {
                size_t s = first, e = last;

                while ( s < e ) {
                    algorithm::swap( array[ s ], array[ e ] );
                    ++s;
                    --e;
                }
            }
        }

        /**
        * @brief Rotates @a k elements of the array. Rotation direction
        * is determined by the sign of @a k, the means a positive value performs
        * a right-rotation and a negative value a left-rotation.
        * @param[in,out] array The array to rotate.
        * @param[in] k Positions to rotate. Sign determines the rotate direction.
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
         * @brief Assigns @a value to all the elements of the array in the
         * range [first,last).
         * @param[in,out] array The array to fill.
         * @param[in] value The value to set all elements to.
         * @param[in] first Initial position of the portion to fill
         * @param[in] last Final position of the portion to fill
         */
        template<typename T, size_t n>
        inline void fill( T ( &array )[ n ],
                         const T value,
                         const size_t first = 0U,
                         const size_t last = n - 1U ) noexcept
        {
            for ( size_t i = first ; i <= last; ++i ) {
                array[ i ] = value;
            }
        }

        /**
        * @brief Returns a pointer to the first element in the range [first,last)
        * that compares equal to @a key. If no such element is found, the
        * function returns @c nullptr.
        * @note The elements are compared using operator '=='
        * @param[in] array The array where the search is performed
        * @param[in] key Value to search for in the range. T shall be a type
        * supporting comparisons using operator==.
        * @param[in] first Initial position of the portion to search
        * @param[in] last Final position of the portion to search
        * @return This function returns a pointer to an entry in the array that
        * matches the search key. If key is not found, a @c nullptr pointer is
        * returned.
        */
        template<typename T, size_t n>
        inline T* find( T ( &array )[ n ],
                        const T key,
                        const size_t first = 0U,
                        const size_t last = n - 1U ) noexcept
        {
            T* found = nullptr;

            for ( size_t i = first; i <= last; ++i ) {
                if ( array[ i ] == key ) {
                    found = &array[ i ];
                    break;
                }
            }
            return found;
        }

        /**
        * @brief Returns @c true if @a pred returns @c true for any of the
        * elements in the range [first,last), and @c false otherwise.
        * @param[in] array The array where the check is performed
        * @param[in] pred Unary function that accepts an element in the range as
        * argument and returns a value convertible to bool. The value returned
        * indicates whether the element fulfills the condition checked by this
        * function.
        * @param[in] first Initial position of the portion to check
        * @param[in] last Final position of the portion to check
        * @return @c true if @a pred returns true for any of the elements in
        * the range [first,last), and @c false otherwise.
        */
        template<typename T, size_t n>
        inline bool any_of( T ( &array )[ n ],
                            bool (*pred)( const T ),
                            const size_t first = 0U,
                            const size_t last = n - 1U ) noexcept
        {
            bool ret = false;

            for ( size_t i = first; i <= last; ++i ) {
                if ( pred( array[ i ] ) ) {
                    ret = true;
                    break;
                }
            }
            return ret;
        }

        /**
        * @brief Returns @c true if @a pred returns @c true for all the
        * elements in the range [first,last), and @c false otherwise.
        * @param[in] array The array where the check is performed
        * @param[in] pred Unary function that accepts an element in the range as
        * argument and returns a value convertible to bool. The value returned
        * indicates whether the element fulfills the condition checked by this
        * function.
        * @param[in] first Initial position of the portion to check
        * @param[in] last Final position of the portion to check
        * @return @c true if @a pred returns true for all the elements in
        * the range [first,last), and @c false otherwise.
        */
        template<typename T, size_t n>
        inline bool all_of( T ( &array )[ n ],
                            bool (*pred)( const T ),
                            const size_t first = 0U,
                            const size_t last = n - 1U ) noexcept
        {
            bool ret = true;

            for ( size_t i = first; i <= last; ++i ) {
                if ( !pred( array[ i ] ) ) {
                    ret = false;
                    break;
                }
            }
            return ret;
        }

        /**
        * @brief Returns the number of elements in the range [first,last) for
        * which @c pred is @c true.
        * @param[in] array The array where the count will be performed
        * @param[in] pred Unary function that accepts an element in the range as
        * argument, and returns a value convertible to bool. The value returned
        * indicates whether the element is counted by this function.
        * @param[in] first Initial position of the portion to check
        * @param[in] last Final position of the portion to check
        * @return The number of elements in the range [first,last) for which
        * @a pred does not return @c false.
        */
        template<typename T, size_t n>
        inline size_t count_if( T ( &array )[ n ],
                                bool (*pred)( const T ),
                                const size_t first = 0U,
                                const size_t last = n - 1U ) noexcept
        {
            size_t count = 0U;

            for ( size_t i = first; i <= last; ++i ) {
                if ( pred( array[ i ] ) ) {
                    ++count;
                }
            }
            return count;
        }

        /**
        * @brief Returns an iterator to the first element in the range [first,last)
        * for which @a pred returns @c true. If no such element is found, the
        * function returns @c nullptr.
        * @param[in] array The array where the search is performed
        * @param[in] pred Unary function that accepts an element in the range as
        * argument and returns a value convertible to bool. The value returned
        * indicates whether the element is considered a match in the context of
        * this function.
        * @param[in] first Initial position of the portion to check
        * @param[in] last Final position of the portion to check
        * @return A pointer to the first element in the range for which @a pred
        * does not return @c false. If @a pred is @c false for all elements,
        * the function returns @c nullptr.
        */
        template<typename T, size_t n>
        inline T* find_if( T ( &array )[ n ],
                           bool (*pred)( const T ),
                           const size_t first = 0U,
                           const size_t last = n - 1U ) noexcept
        {
            T *found = nullptr;

            for ( size_t i = first; i <= last; ++i ) {
                if ( pred( array[ i ] ) ) {
                    found = &array[ i ];
                }
            }
            return found;
        }

        /**
        * @brief Applies function @a fn to each of the elements in the range
        *  [first,last).
        * @param[in] array The array
        * @param[in] fn Unary function that accepts an element in the range as
        *  argument.
        * @param[in] first Initial position of the portion to check
        * @param[in] last Final position of the portion to check
        * @return none
        */
        template<typename T, size_t n>
        inline void for_each( T ( &array )[ n ],
                              void (*fn)( T& ),
                              const size_t first = 0U,
                              const size_t last = n - 1U ) noexcept
        {
            for ( size_t i = first; i <= last; ++i ) {
                (void)fn( array[ i ] );
            }
        }

        /**
        * @brief Returns a pointer to the first element in the range [first,last)
        * that compares equal to @a key. If no such element is found, the
        * function returns @c nullptr.
        * @note The elements in the range shall already be sorted according to
        * this same criterion (operator< or operator==), or at least partitioned with
        * respect to @a key.
        * @note The elements are compared using operator '=='
        * @param[in] array The array where the search is performed
        * @param[in] key Value to search for in the range. T shall be a type
        * supporting comparisons using operator== and operator<.
        * @param[in] first Initial position of the portion to search
        * @param[in] last Final position of the portion to search
        * @return This function returns a pointer to an entry in the array that
        * matches the search key. If key is not found, a @c nullptr pointer is
        * returned.
        */
        template<typename T, size_t n>
        inline T* binary_search( T ( &array )[ n ],
                                 const T key,
                                 const size_t first = 0U,
                                 const size_t last = n - 1U ) noexcept
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