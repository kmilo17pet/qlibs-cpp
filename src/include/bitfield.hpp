#ifndef QLIBS_BITFIELD
#define QLIBS_BITFIELD

#include "include/types.hpp"

namespace qlibs {

    constexpr size_t bitfield_size( const size_t nbits )
    {
        return 4U*( ( ( nbits - 1U )/32U ) + 1U );
    }

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
            bool setup( void * const area,
                        const size_t area_size ) noexcept;
            bool clearAll( void ) noexcept;
            bool setAll( void ) noexcept;
            bool setBit( const size_t index ) noexcept;
            bool clearBit( const size_t index ) noexcept;
            bool toggleBit( const size_t index ) noexcept;
            bool readBit( const size_t index ) const noexcept;
            bool writeBit( const size_t index,
                           const bool value ) noexcept;
            uint32_t readUINTn( const size_t index,
                                const size_t xBits ) const noexcept;
            bool writeUINTn( const size_t index,
                             const size_t xBits,
                             uint32_t value ) noexcept;
            float readFloat( const size_t index ) const noexcept;
            bool writeFloat( const size_t index,
                             const float value ) noexcept;
            void* dump( void * const dst,
                        const size_t n ) noexcept;
    };

}

#endif /*QLIBS_BITFIELD*/