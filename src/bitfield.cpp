#include <include/bitfield.hpp>

using namespace qlibs;

const size_t bitfield::LBit =  static_cast<size_t>( sizeof(uint32_t) * 8U );

/*============================================================================*/
bool bitfield::setup( void * const area,
                      const size_t area_size ) noexcept
{
    bool retValue = false;

    if ( ( nullptr != area ) && ( area_size > 0U ) ) {
        /*cstat -CERT-EXP36-C_b*/
        field = static_cast<uint32_t *>( area );
        /*cstat +CERT-EXP36-C_b*/
        size = area_size*8U;
        nSlots = area_size/sizeof(uint32_t);
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool bitfield::clearAll( void ) noexcept
{
    bool retValue = false;

    if ( nullptr != field ) {
        (void)memset( field, 0, size/8U );
        retValue = true;
    }
    return retValue;
}
/*============================================================================*/
bool bitfield::setAll( void ) noexcept
{
    bool retValue = false;

    if ( nullptr != field ) {
        (void)memset( field, 0xFF, size/8U );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool bitfield::setBit( const size_t index ) noexcept
{
    bool retValue = false;

    if ( nullptr != field ) {
        set( index );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool bitfield::clearBit( const size_t index ) noexcept
{
    bool retValue = false;

    if ( nullptr != field ) {
        clear( index );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool bitfield::toggleBit( const size_t index ) noexcept
{
    bool retValue = false;

    if ( nullptr != field ) {
        toggle( index );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
bool bitfield::readBit( const size_t index ) const noexcept
{
    bool retValue = false;

    if ( nullptr != field ) {
        retValue = 0U != get( index );
    }

    return retValue;
}
/*============================================================================*/
bool bitfield::writeBit( const size_t index,
                         const bool value ) noexcept
{
    bool retValue = false;

    if ( nullptr != field ) {
        if ( value ) {
            set( index );
        }
        else {
            clear( index );
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
uint32_t bitfield::readUINTn( const size_t index,
                              const size_t xBits ) const noexcept
{
    uint32_t retValue = 0U;

    if ( ( nullptr != field ) && ( xBits <= 32U ) ) {
        if ( 1U == xBits ) {
            retValue = static_cast<uint32_t>( readBit( index ) );
        }
        else if ( 32U == xBits ) {
            retValue = read_uint32( index );
        }
        else {
            retValue = read_uint32( index );
            retValue &= safeMask( 0xFFFFFFFFU, 32U, xBits );
        }
    }

    return retValue;
}
/*============================================================================*/
bool bitfield::writeUINTn( const size_t index,
                           const size_t xBits,
                           uint32_t value ) noexcept
{
    bool retValue = false;

    if ( ( nullptr != field ) && ( xBits <= 32U ) ) {
        uint32_t w, wMask;

        if ( 1U == xBits ) {
            (void)writeBit( index, 0 != value );
        }
        else if ( 32U == xBits ) {
            write_uint32( index, value );
        }
        else {
            w = read_uint32( index );
            value &= safeMask( 0xFFFFFFFFU, 32U, xBits );
            /*cstat -ATH-overflow*/
            wMask = static_cast<uint32_t>( 0xFFFFFFFFU ) << static_cast<uint32_t>( xBits );
            /*cstat +ATH-overflow*/
            write_uint32( index, maskMerge( w, value, wMask ) );
        }
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
float bitfield::readFloat( const size_t index ) const noexcept
{
    float retValue = 0.0F;

    if ( nullptr != field ) {
        uint32_t rval;

        rval = read_uint32( index );
        (void)memcpy( &retValue, &rval, sizeof(float) );
    }

    return retValue;
}
/*============================================================================*/
bool bitfield::writeFloat( const size_t index,
                           const float value ) noexcept
{
    bool retValue = false;

    if ( nullptr != field ) {
        uint32_t fval = 0U;

        (void)memcpy( &fval, &value, sizeof(float) );
        write_uint32( index, fval );
        retValue = true;
    }

    return retValue;
}
/*============================================================================*/
void* bitfield::dump( void * const dst,
                      const size_t n ) noexcept
{
    void *retValue = nullptr;

    if ( ( nullptr != field ) && ( nullptr != dst ) ) {
        if ( n <= ( size/8U ) ) {
            retValue = memcpy( dst, static_cast<const void*>( field ), n );
        }
    }

    return retValue;
}
/*============================================================================*/
uint32_t bitfield::read_uint32( const size_t index ) const noexcept
{
    size_t s, of, bits_taken;
    uint32_t result;

    s = slot( index );
    of = offset( index );
    /*cstat -CERT-INT34-C_a*/
    result = field[ s ] >> of;
    /*cstat +CERT-INT34-C_a*/
    bits_taken  = 32U - of;
    if ( ( 0U != of ) && ( ( index + bits_taken ) < size ) ) {
        /*cstat -ATH-shift-bounds -MISRAC++2008-5-8-1 -CERT-INT34-C_b*/
        result |= field[ s + 1U ] << static_cast<uint32_t>( bits_taken );
        /*cstat +ATH-shift-bounds +MISRAC++2008-5-8-1 +CERT-INT34-C_b*/
    }

    return result;
}
/*============================================================================*/
void bitfield::write_uint32( const size_t index,
                             const uint32_t value ) noexcept
{
    uint32_t wMask;
    size_t s, of;

    s = slot( index );
    of = offset( index );
    if ( 0U == of ) {
        field[ s ] = value;
    }
    else {
        wMask = safeMask( 0xFFFFFFFFU, LBit, of );
        /*cstat -CERT-INT34-C_a*/
        field[ s ] = ( value << of ) | ( field[ s ] & wMask );
        /*cstat +CERT-INT34-C_a*/
        if ( ++s < nSlots ) {
            field[ s ] = safeMask( value, 32U, of ) | ( field[ s ] & ( ~wMask ) );
        }
    }
}
/*============================================================================*/