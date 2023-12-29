/*!
 * @file bitfield.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief A bit-field manipulation library.
 **/


#ifndef QLIBS_BITFIELD
#define QLIBS_BITFIELD

#include <include/qlibs_types.hpp>

/**
* @brief The qLibs++ library namespace.
*/
namespace qlibs {

    /** @addtogroup qbitfield Bit-Field
    * @brief API for the Bit-Field manipulation library
    *  @{
    */

    /**
    * @brief Determine the @c uint8_t array-size to setup a BitField instance.
    * @param[in] n The desired number of bits for the BitField.
    * @return The number for bytes ( or array size for @c uint8_t ) required for
    * the BitField of @a n bits
    */
    constexpr size_t bitfield_size( const size_t n )
    {
        return 4U*( ( ( n - 1U )/32U ) + 1U );
    }

    /**
    * @brief A BitField object
    */
    class bitfield : private nonCopyable {
        private:
            uint32_t *field{ nullptr };
            size_t size{ 0U };
            size_t nSlots{ 0U };
            static const size_t LBit;
            inline uint32_t mask( const size_t index ) noexcept
            {
                return static_cast<uint32_t>( 1U ) << ( index % LBit );
            }
            inline size_t bitSlot( const size_t index ) const noexcept
            {
                return index/LBit;
            }
            inline uint32_t bitGet( const size_t index ) const noexcept
            {
                const size_t slot = bitSlot( index );
                return  ( field[ slot ] >> ( index % LBit ) ) & 1U;
            }
            inline void bitSet( const size_t index ) noexcept
            {
                const size_t slot = bitSlot( index );
                field[ slot ] |= mask( index );
            }
            inline void bitClear( const size_t index ) noexcept
            {
                const size_t slot = bitSlot( index );

                field[ slot ] &= ~mask( index );
            }
            inline void bitToggle( const size_t index ) noexcept
            {
                const size_t slot = bitSlot( index );

                field[ slot ] ^= mask( index );
            }
            inline uint32_t safeMask( const uint32_t val,
                                      const size_t x,
                                      const size_t nbits ) const noexcept
            {
                return val >> ( static_cast<uint32_t>( x - nbits ) );
            }
            inline size_t offset( const size_t index ) const noexcept
            {
                return index & static_cast<size_t>( 31U );
            }
            inline uint32_t maskMerge( const uint32_t w,
                                       const uint32_t value,
                                       const uint32_t mask ) noexcept
            {
                return value ^ ( ( w ^ value ) & mask );
            }

            uint32_t read_uint32( const size_t index ) const noexcept;
            void write_uint32( const size_t index,
                               const uint32_t value ) noexcept;

        public:
            bitfield() = default;
            virtual ~bitfield() {}

            /**
            * @brief Setup a initialize a BitField instance.
            * @param[in] area A pointer to the memory block to hold the BitField.
            * Should be an uint8_t array of size bitfield_size(n), where n, is the
            * number of bits inside the BitField.
            * @param[in] area_size The number of bytes in @a area.
            * @return @c true on success, otherwise return @c false.
            */
            bool setup( void * const area,
                        const size_t area_size ) noexcept;

            /**
            * @brief Clear all the bits in the BitField.
            * @return @c true on success, otherwise return @c false.
            */
            bool clearAll( void ) noexcept;

            /**
            * @brief Set all the bits in the BitField.
            * @return @c true on success, otherwise return @c false.
            */
            bool setAll( void ) noexcept;

            /**
            * @brief Sets one bit in a BitField
            * @param[in] index The bit-index.
            * @return @c true on success, otherwise return @c false.
            */
            bool setBit( const size_t index ) noexcept;

            /**
            * @brief Clears one bit in a BitField
            * @param[in] index The bit-index.
            * @return @c true on success, otherwise return @c false.
            */
            bool clearBit( const size_t index ) noexcept;

            /**
            * @brief Toggles (i.e. reverses the state of) a bit in a BitField
            * @param[in] index The bit-index.
            * @return @c true on success, otherwise return @c false.
            */
            bool toggleBit( const size_t index ) noexcept;

            /**
            * @brief Retrieve the state of a bit in a bitfield
            * @param[in] index The bit-index.
            * @return The value of the bit at @a index.
            */
            bool readBit( const size_t index ) const noexcept;

            /**
            * @brief Writes one bit in a bitfield
            * @param[in] index The bit-index.
            * @param[in] value The boolean value to write.
            * @return @c true on success, otherwise return @c false.
            */
            bool writeBit( const size_t index,
                           const bool value ) noexcept;

            /**
            * @brief Reads an unsigned 32-bit value from the BitField
            * @param[in] index The bit-index taken as offset.
            * @param[in] xBits The number of bits to read. ( max allowed : 32 bits )
            * @return The value from the bitfield from the desired index
            */
            uint32_t readUINTn( const size_t index,
                                const size_t xBits ) const noexcept;

            /**
            * @brief Writes an unsigned n-bit value from the BitField
            * @param[in] index The bit-index taken as offset.
            * @param[in] value The value to write.
            * @param[in] xBits The number of bits to read. ( max allowed : 32 bits )
            * @return @c true on success, otherwise return @c false.
            */
            bool writeUINTn( const size_t index,
                             const size_t xBits,
                             uint32_t value ) noexcept;

            /**
            * @brief Reads a 32-bit floating point value from the BitField
            * @param[in] index The bit-index taken as offset.
            * @return The floating point value from the BitField at the desired index
            */
            float readFloat( const size_t index ) const noexcept;

            /**
            * @brief Writes a 32-bit floating point value to the BitField
            * @param[in] index The bit-index taken as offset.
            * @param[in] value The floating point value to write.
            * @return @c true on success, otherwise return @c false.
            */
            bool writeFloat( const size_t index,
                             const float value ) noexcept;

            /**
            * @brief Copies @a n bytes from the bit-field instance to a designed memory
            * area.
            * @param[in] dst Pointer to the destination array where the content is to
            * be copied.
            * @param[in] n Number of bytes to copy.
            * @return Destination is returned on success, otherwise @c nullptr.
            */
            void* dump( void * const dst,
                        const size_t n ) noexcept;
    };

    /** @}*/
}


#endif /*QLIBS_BITFIELD*/