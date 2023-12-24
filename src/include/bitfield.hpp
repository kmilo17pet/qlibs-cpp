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
            inline uint32_t mask( const size_t index )
            {
                return static_cast<uint32_t>( 1U ) << ( index % LBit );
            }
            inline size_t bitSlot( const size_t index ) const
            {
                return index/LBit;
            }
            inline uint32_t bitGet( const size_t index ) const
            {
                const size_t slot = bitSlot( index );
                return  ( field[ slot ] >> ( index % LBit ) ) & 1U;
            }
            inline void bitSet( const size_t index )
            {
                const size_t slot = bitSlot( index );
                field[ slot ] |= mask( index );
            }
            inline void bitClear( const size_t index )
            {
                const size_t slot = bitSlot( index );

                field[ slot ] &= ~mask( index );
            }
            inline void bitToggle( const size_t index )
            {
                const size_t slot = bitSlot( index );

                field[ slot ] ^= mask( index );
            }
            inline uint32_t safeMask( const uint32_t val, const size_t x, const size_t nbits ) const
            {
                return val >> ( static_cast<uint32_t>( x - nbits ) );
            }
            inline size_t offset( const size_t index ) const
            {
                return index & static_cast<size_t>( 31U );
            }
            inline uint32_t maskMerge( const uint32_t w, const uint32_t value,  const uint32_t mask )
            {
                return value ^ ( ( w ^ value ) & mask );
            }

            uint32_t read_uint32( const size_t index ) const;
            void write_uint32( const size_t index, const uint32_t value );

        public:
            bitfield() = default;
            virtual ~bitfield() {}
            bool setup( void * const area, const size_t area_size );
            bool clearAll( void );
            bool setAll( void );
            bool setBit( const size_t index );
            bool clearBit( const size_t index );
            bool toggleBit( const size_t index );
            bool readBit( const size_t index ) const;
            bool writeBit( const size_t index, bool value );
            uint32_t readUINTn( const size_t index, size_t xBits ) const;
            bool writeUINTn( const size_t index, size_t xBits, uint32_t value );
            float readFloat( const size_t index ) const;
            bool writeFloat( const size_t index, float value );
            void* dump( void * const dst, size_t n );
    };

}

#endif /*QLIBS_BITFIELD*/