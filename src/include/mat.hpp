/*!
 * @file mat.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief Computes the RMS (Root Mean Square) of a signal using a 2-step
 * recursive average specially designed for micro-controllers with FPU.
 **/

#ifndef QLIBS_MAT
#define QLIBS_MAT

#include <include/qlibs_types.hpp>
#include <include/ffmath.hpp>

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {

    /** @addtogroup qmat Matrix library
    * @brief Class for Recursive Root Mean Square RMS estimation
    *  @{
    */

    enum matSpecial {
        MAT_ZERO,
        MAT_IDENTITY,
        MAT_ALL,
        MAT_MAGIC,
    };

    /**
    * @brief Matrix object
    */
// Fixed-size matrix class without dynamic memory allocation
    template <size_t Rows, size_t Cols>
    class mat {
        template <size_t R, size_t C>
        friend class mat;
        private:
            real_t m[Rows][Cols]; // Statically allocated memory for the matrix
        public:
            mat()
            {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        m[i][j] = 0.0_re; // Initialize matrix elements to zero
                    }
                }
            }

            mat( const matSpecial type, const real_t value = 1.0_re )
            {
                switch( type ) {
                    case MAT_IDENTITY:
                        static_assert(Rows == Cols, "Identity matrix must be square");

                        for (size_t i = 0; i < Rows; ++i) {
                            for (size_t j = 0; j < Cols; ++j) {
                                if (i == j) {
                                    m[i][j] = value;
                                } else {
                                    m[i][j] = 0.0_re;
                                }
                            }
                        }
                        break;
                    case MAT_ALL:
                        for (size_t i = 0; i < Rows; ++i) {
                            for (size_t j = 0; j < Cols; ++j) {
                                m[i][j] = value;
                            }
                        }
                        break;
                    case MAT_MAGIC:

                        break;
                    default:
                        for (size_t i = 0; i < Rows; ++i) {
                            for (size_t j = 0; j < Cols; ++j) {
                                m[i][j] = 0.0_re; // Initialize matrix elements to zero
                            }
                        }
                        break;
                }
            }

            template <typename... Args>
            mat(Args... args)
            {
                static_assert(sizeof...(args) == Rows * Cols, "Incorrect number of elements provided");
                real_t elements[] = {static_cast<real_t>(args)...};
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        m[i][j] = elements[i * Cols + j];
                    }
                }
            }

            //inline real_t& operator()(int row, int col)
            //{
            //    return m[row][col];
            //}

            inline real_t& operator()(size_t row, size_t col)
            {
                return m[row][col];
            }

            inline real_t& operator()(size_t index)
            {
                return m[index / Cols][index % Cols];
            }

            inline const real_t& operator()(size_t index) const
            {
                return m[index / Cols][index % Cols];
            }

            // Copy constructor declaration and definition
            mat(const mat<Rows, Cols>& other)
            {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        m[i][j] = other.m[i][j]; // Copy elements
                    }
                }
            }

            // Overload the + operator to perform matrix addition
            mat<Rows, Cols> operator+(const mat<Rows, Cols>& other) const
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = m[i][j] + other.m[i][j];
                    }
                }
                return result;
            }
            // Overload the subtraction operator to perform matrix subtraction
            mat<Rows, Cols> operator-(const mat<Rows, Cols>& other) const
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = m[i][j] - other.m[i][j];
                    }
                }
                return result;
            }
            // Overload the * operator to perform matrix multiplication
            template <size_t N>
            mat<Rows, N> operator*(const mat<Cols, N>& other) const
            {
                mat<Rows, N> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < N; ++j) {
                        result.m[i][j] = 0;
                        for (size_t k = 0; k < Cols; ++k) {
                            result.m[i][j] += m[i][k] * other.m[k][j];
                        }
                    }
                }
                return result;
            }

            // Overload the multiplication operator for scalar multiplication
            mat<Rows, Cols> operator*(const real_t& scalar) const
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = m[i][j] * scalar;
                    }
                }
                return result;
            }

            // Overload the division operator for scalar multiplication
            mat<Rows, Cols> operator/(const real_t& scalar) const
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = m[i][j] / scalar;
                    }
                }
                return result;
            }

            // Overload the multiplication operator for scalar * matrix
            friend mat<Rows, Cols> operator*(const real_t& scalar, const mat<Rows, Cols>& matrix)
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = scalar * matrix.m[i][j];
                    }
                }
                return result;
            }

            // Overload the divison operator for scalar / matrix
            //friend mat<Rows, Cols> operator/(const real_t& scalar, const mat<Rows, Cols>& matrix)
            //{
            //    mat<Rows, Cols> result;
            //    for (size_t i = 0; i < Rows; ++i) {
            //        for (size_t j = 0; j < Cols; ++j) {
            //            result.m[i][j] = scalar / matrix.m[i][j];
            //        }
            //    }
            //    return result;
            //}

            mat<Rows, Cols> operator+(const real_t& scalar) const
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = m[i][j] + scalar;
                    }
                }
                return result;
            }

            // Overload the unary increment operator
            mat<Rows, Cols>& operator++()
            {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        ++m[i][j];
                    }
                }
                return *this;
            }
            // Overload the post-increment operator (++)
            mat<Rows, Cols> operator++(int)
            {
                mat<Rows, Cols> temp = *this; // Make a copy of the current matrix
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        ++m[i][j];
                    }
                }
                return temp; // Return the matrix before incrementing
            }

            mat<Rows, Cols> operator-(const real_t& scalar) const
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = m[i][j] - scalar;
                    }
                }
                return result;
            }

            // Overload the unary increment operator
            mat<Rows, Cols>& operator--()
            {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        ++m[i][j];
                    }
                }
                return *this;
            }
            // Overload the post-increment operator (++)
            mat<Rows, Cols> operator--(int)
            {
                mat<Rows, Cols> temp = *this; // Make a copy of the current matrix
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        --m[i][j];
                    }
                }
                return temp; // Return the matrix before incrementing
            }

            // Overload the compound addition operator (+=) for matrices
            mat<Rows, Cols>& operator+=(const mat<Rows, Cols>& other) {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        m[i][j] += other.m[i][j];
                    }
                }
                return *this;
            }

            // Overload the compound subtraction operator (-=) for matrices
            mat<Rows, Cols>& operator-=(const mat<Rows, Cols>& other) {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        m[i][j] -= other.m[i][j];
                    }
                }
                return *this;
            }

            // Overload the compound multiplication operator (*=) with another matrix
            template <size_t N>
            mat<Rows, N>& operator*=(const mat<Cols, N>& other) {
                mat<Rows, N> result;

                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < N; ++j) {
                        result.m[i][j] = 0;
                        for (size_t k = 0; k < Cols; ++k) {
                            result.m[i][j] += m[i][k] * other.m[k][j];
                        }
                    }
                }

                *this = result; // Update the current matrix with the result
                return *this;
            }

            // Overload the compound multiplication operator (*=) with scalar
            mat<Rows, Cols>& operator*=(const real_t& scalar) {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        m[i][j] *= scalar;
                    }
                }
                return *this;
            }

            // Overload the compound multiplication operator (*=) with scalar
            mat<Rows, Cols>& operator/=(const real_t& scalar) {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        m[i][j] /= scalar;
                    }
                }
                return *this;
            }

            // Overload the compound addition operator (+=) with a scalar
            mat<Rows, Cols>& operator+=(const real_t& scalar)
            {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        m[i][j] += scalar;
                    }
                }
                return *this;
            }

            // Overload the compound subtraction operator (-=) with a scalar
            mat<Rows, Cols>& operator-=(const real_t& scalar)
            {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        m[i][j] -= scalar;
                    }
                }
                return *this;
            }

            // Overload the addition operator for scalar + matrix
            friend mat<Rows, Cols> operator+(const real_t& scalar, const mat<Rows, Cols>& matrix)
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = scalar + matrix.m[i][j];
                    }
                }
                return result;
            }

            // Overload the addition operator for scalar + matrix
            friend mat<Rows, Cols> operator-(const real_t& scalar, const mat<Rows, Cols>& matrix)
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = scalar - matrix.m[i][j];
                    }
                }
                return result;
            }

            // Overload the assignment operator
            mat<Rows, Cols>& operator=(const mat<Rows, Cols>& other)
            {
                if (this != &other) { // Avoid self-assignment
                    for (size_t i = 0; i < Rows; ++i) {
                        for (size_t j = 0; j < Cols; ++j) {
                            m[i][j] = other.m[i][j]; // Copy elements
                        }
                    }
                }
                return *this;
            }

            // Overload the unary minus operator
            mat<Rows, Cols> operator-() const
            {
                mat<Rows, Cols> result;
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        result.m[i][j] = -m[i][j]; // Negate elements
                    }
                }
                return result;
            }

            // Display matrix elements
            void display() const
            {
                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        std::cout << m[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            }

            mat<Cols, Rows> operator!() const {
                mat<Cols, Rows> transposed;

                for (size_t i = 0; i < Rows; ++i) {
                    for (size_t j = 0; j < Cols; ++j) {
                        transposed[j][i] = m[i][j];
                    }
                }

                return transposed;
            }

            // Overload the [][] operator to access elements of the matrix
            real_t* operator[](size_t index)
            {
                return m[index];
            }

            const real_t* operator[](size_t index) const
            {
                return m[index];
            }

            mat<Rows, Rows> inv() const
            {
                static_assert(Rows == Cols, "Inverse exists only for square matrices");
                mat<Rows,Rows> X = *this;
                mat<Rows,Rows> I( MAT_IDENTITY );

                for ( size_t i = 0; i < Rows; ++i ) {
                    real_t pivot = X(i, i);
                    // Make the diagonal element 1
                    for ( size_t j = 0; j < Rows; ++j ) {
                        X(i, j) /= pivot;
                        I(i, j) /= pivot;
                    }
                    for ( size_t j = 0; j < Rows; ++j) { // Make other elements in the column 0
                        if ( j != i ) {
                            real_t factor = X(j, i);
                            for ( size_t k = 0; k < Rows; ++k ) {
                                X(j, k) -= factor * X(i, k);
                                I(j, k) -= factor * I(i, k);
                            }
                        }
                    }
                }

                return I;
            }

            real_t normInf() const
            {
                size_t i,j;
                real_t maxM =0.0_re, sum = 0.0_re;

                for ( i = 0U ; i < Cols; i++ ) {
                    sum = 0.0_re;
                    for ( j = 0U ; j < Rows ; ++j ) {
                        sum += ( m[i][j] < 0.0_re ) ? -m[i][j] : m[i][j];
                    }
                    if ( 0U == i ) {
                        maxM = sum;
                    }
                    maxM = ( sum > maxM ) ? sum : maxM;
                }
                return maxM;
            }
    };



    /** @}*/
}


#endif /*QLIBS_MAT*/