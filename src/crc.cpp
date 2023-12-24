#include <include/crc.hpp>

using namespace qlibs;

/*============================================================================*/
uint32_t crc::reflect( uint32_t xData, const uint8_t nBits )
{
    uint32_t  r = 0;
    uint8_t xBit;
    /*Reflect the data about the center bit*/
    for ( xBit= 0U ; xBit < nBits ; ++xBit ) {
        /*if the LSB bit is set, set the reflection of it*/
        if ( 0U != ( xData & 0x01U ) ) {
            /*cstat -MISRAC++2008-5-0-10 -MISRAC++2008-5-0-8*/
            r |= static_cast<uint32_t>( 1U << ( ( nBits - 1U ) - xBit ) );
            /*cstat +MISRAC++2008-5-0-10 +MISRAC++2008-5-0-8*/
        }
        xData >>= 1U;
    }

    return r;
}
/*============================================================================*/
uint32_t crc::generic( crcMode mode, const void * const pData, const size_t length, uint32_t poly, const uint32_t init, bool refIn, bool refOut, uint32_t xorOut )
{
    uint32_t val = 0U;

    if ( ( nullptr != pData ) && ( length > 0U ) && ( ( mode == crcMode::CRC8 ) || ( mode == crcMode::CRC16 ) || ( mode == crcMode::CRC32 ) ) ) {
        size_t i;
        uint8_t xBit;
        const uint32_t widthValues[ 3 ] = { 8UL, 16UL, 32UL };
        const uint32_t width = widthValues[ mode ];
        /*cstat -CERT-EXP36-C_b*/
        uint8_t const * const msg = static_cast<uint8_t const *>( pData );
        /*cstat +CERT-EXP36-C_b*/
        const uint32_t wd1 = static_cast<uint32_t>( width ) - 8U;
        const uint32_t topBit = static_cast<uint32_t>( 1U ) << ( width - 1U );
        const uint32_t bitMask = ( 0xFFFFFFFFU >> ( 32U - width ) );
        poly &= bitMask;
        xorOut &= bitMask;
        val = init;

        for ( i = 0 ; i < length ; ++i ) {
            /*cstat -MISRAC++2008-5-0-8*/
            val ^= ( refIn ) ? ( reflect( static_cast<uint32_t>( msg[ i ] ) , 8U ) <<  wd1 )
                             : ( static_cast<uint32_t>( msg[ i ] ) << wd1 );
            /*cstat +MISRAC++2008-5-0-8*/
            for ( xBit = 8U ; xBit > 0U ; --xBit ) {
                val = ( 0U != ( val & topBit ) ) ? ( ( val << 1U ) ^ poly )
                                                 : ( val << 1U );
            }
        }
        val = ( refOut ) ? ( reflect( val, static_cast<uint8_t>( width ) )^xorOut ) : ( val^xorOut );
        val &= bitMask;
    }

    return val;
}
/*============================================================================*/