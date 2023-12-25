#ifndef QLIBS_CRC
#define QLIBS_CRC

#include "include/types.hpp" 


namespace qlibs {

    enum crcMode {
        CRC8 = 0,
        CRC16,
        CRC32,
    };

    class crc {
        private:
            static uint32_t reflect( uint32_t xData,
                                     const uint8_t nBits ) noexcept;
        public:
            virtual ~crc() {}
            crc() = default;
            uint32_t generic( crcMode mode,
                              const void * const pData,
                              const size_t length,
                              uint32_t poly,
                              const uint32_t init = 0U,
                              bool refIn = false,
                              bool refOut = false,
                              uint32_t xorOut = 0U ) noexcept;
            inline uint8_t crc8( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x07U ) );
            }
            inline uint8_t crc8_CDMA2000( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x9BU, 0xFFU ) );
            }
            inline uint8_t crc8_DARC( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x39U, 0U, true, true ) );
            }
            inline uint8_t crc8_DVS_S2( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0xD5U ) );
            }
            inline uint8_t crc8_EBU( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x1DU, 0xFFU, true, true ) );
            }
            inline uint8_t crc8_I_CODE( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x1DU, 0xFDU ) );
            }
            inline uint8_t crc8_ITU( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x07U, 0U, false, false, 0x55U ) );
            }
            inline uint8_t crc8_MAXIM( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x31U, 0U, true, true ) );
            }
            inline uint8_t crc8_ROHC( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x07U, 0xFFU, true, true ) );
            }
            inline uint8_t crc8_WCDMA( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x9BU, 0U, true, true ) );
            }
            inline uint16_t crc16_CCITT_FALSE( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xFFFFU ) );
            }
            inline uint16_t crc16_ARC( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0U, true, true ) );
            }
            inline uint16_t crc16_AUG_CCITT( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0x1D0FU ) );
            }
            inline uint16_t crc16_BUYPASS( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U ) );
            }
            inline uint16_t crc16_CDMA2000( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0xC867U, 0xFFFFU ) );
            }
            inline uint16_t crc16_DDS_110( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0x800DU ) );
            }
            inline uint16_t crc16_DECT_R( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x0589U, 0U, false, false, 0x0001U ) );
            }
            inline uint16_t crc16_DECT_X( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x0589U ) );
            }
            inline uint16_t crc16_DNP( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x3D65U, 0U, true, true, 0xFFFFU ) );
            }
            inline uint16_t crc16_EN_13757( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x3D65U, 0U, false, false, 0xFFFFU ) );
            }
            inline uint16_t crc16_GENIBUS( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xFFFFU, false, false, 0xFFFFU ) );
            }
            inline uint16_t crc16_MAXIM( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0U, true, true, 0xFFFFU ) );
            }
            inline uint16_t crc16_MCRF4XX( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xFFFFU, true, true ) );
            }
            inline uint16_t crc16_RIELLO( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xB2AAU, true, true ) );
            }
            inline uint16_t crc16_T10_DIF( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8BB7U ) );
            }
            inline uint16_t crc16_TELEDISK( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0xA097U ) );
            }
            inline uint16_t crc16_TMS37157( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0x89ECU, true, true ) );
            }
            inline uint16_t crc16_USB( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0xFFFFU, true, true, 0xFFFFU ) );
            }
            inline uint16_t crc16_A( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xc6c6U, true, true ) );
            }
            inline uint16_t crc16_KERMIT( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0U, true, true ) );
            }
            inline uint16_t crc16_MODBUS( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0xFFFFU, true, true ) );
            }
            inline uint16_t crc16_X_25( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xFFFFU, true, true, 0xFFFFU ) );
            }
            inline uint16_t crc16_XMODEM( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U ) );
            }
            inline uint32_t crc32( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0xFFFFFFFFU, true, true, 0xFFFFFFFFU  );
            }
            inline uint32_t crc32_BZIP2( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0xFFFFFFFFU, false, false, 0xFFFFFFFFU  );
            }
            inline uint32_t crc32_C( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x1EDC6F41U, 0xFFFFFFFFU, true, true, 0xFFFFFFFFU  );
            }
            inline uint32_t crc32_D( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0xA833982BU, 0xFFFFFFFFU, true, true, 0xFFFFFFFFU  );
            }
            inline uint32_t crc32_JAMCRC( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0xFFFFFFFFU, true, true );
            }
            inline uint32_t crc32_MPEG2( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0xFFFFFFFFU );
            }
            inline uint32_t crc32_POSIX( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0U, false, false, 0xFFFFFFFFU  );
            }
            inline uint32_t crc32_Q( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x814141ABU );
            }
            inline uint32_t crc32_XFER( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x000000AFU );
            }
    };

}

#endif /*QLIBS_CRC*/