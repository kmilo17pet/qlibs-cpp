/*!
 * @file bitfield.hpp
 * @author J. Camilo Gomez C.
 * @version 1.01
 * @note This file is part of the qLibs-cpp distribution.
 * @brief A bit-field manipulation library.
 **/

#ifndef QLIBS_CRC
#define QLIBS_CRC

#include "include/types.hpp" 

namespace qlibs {

    /** @addtogroup qcrc CRC
    * @brief Generic Cyclic Redundancy Check (CRC) calculator class
    *  @{
    */

    /**
    * @brief Enumeration with all the supported cyclic redundancy checks
    */ 
    enum crcMode {
        CRC8 = 0,    /*!< 8-Bit Cyclic Redundancy Check*/
        CRC16,       /*!< 16-Bit Cyclic Redundancy Check*/
        CRC32,       /*!< 32-Bit Cyclic Redundancy Check*/
    };

    class crc {
        private:
            static uint32_t reflect( uint32_t xData,
                                     const uint8_t nBits ) noexcept;
            crc() = default;
        public:

            /**
            * @brief Calculates in one pass the common @a width bit CRC value for a
            * block of data that is passed to the function together with a parameter
            * indicating the @a length.
            * @param[in] mode To select the CRC calculation mode. Only the following
            * values are supported: ::CRC8, ::CRC16 and ::CRC32.
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @param[in] poly CRC polynomial value.
            * @param[in] init CRC initial value.
            * @param[in] refIn If @c true, the input data is reflected before processing.
            * @param[in] refOut If @c true, the CRC result is reflected before output.
            * @param[in] xorOut The final XOR value.
            * @return The CRC value for @a data.
            */
            static uint32_t generic( crcMode mode,
                                     const void * const pData,
                                     const size_t length,
                                     uint32_t poly,
                                     const uint32_t init = 0U,
                                     bool refIn = false,
                                     bool refOut = false,
                                     uint32_t xorOut = 0U ) noexcept;

            /**
            * @brief CRC-8 with poly = 0x07 init = 0x00 refIn = false refOut = false
            * xorOut= 0x00
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
            */
            static inline uint8_t crc8( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x07U ) );
            }

            /**
            * @brief CRC-8/CDMA2000 with poly = 0x9B init = 0xFF refIn = false
            * refOut = false  xorOut= 0x00
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint8_t crc8_CDMA2000( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x9BU, 0xFFU ) );
            }

            /**
            * @brief CRC-8/DARC with poly = 0x39 init = 0x00 refIn = true
            * refOut = true  xorOut= 0x00
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint8_t crc8_DARC( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x39U, 0U, true, true ) );
            }

            /**
            * @brief CRC-8/DVS-S2 with poly = 0xD5 init = 0x00 refIn = false
            * refOut = false xorOut= 0x00 
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint8_t crc8_DVS_S2( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0xD5U ) );
            }

            /**
            * @brief CRC-8/EBU with poly = 0x1D init = 0xFF refIn = true refOut = true
            * xorOut= 0x00 
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint8_t crc8_EBU( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x1DU, 0xFFU, true, true ) );
            }

            /**
            * @brief CRC-8/I-CODE with poly = 0x1D init = 0xFD refIn = false
            * refOut = false  xorOut= 0x00
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint8_t crc8_I_CODE( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x1DU, 0xFDU ) );
            }

            /**
            * @brief CRC-8/ITU with poly = 0x07 init = 0x00 refIn = false
            * refOut = false  xorOut= 0x55
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint8_t crc8_ITU( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x07U, 0U, false, false, 0x55U ) );
            }

            /**
            * @brief CRC-8/MAXIM with poly = 0x31 init = 0x00 refIn = true 
            * refOut = true  xorOut= 0x00
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint8_t crc8_MAXIM( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x31U, 0U, true, true ) );
            }

            /**
            * @brief CRC-8/ROHC with poly = 0x07 init = 0xFF refIn = true 
            * refOut = true  xorOut= 0x00
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint8_t crc8_ROHC( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x07U, 0xFFU, true, true ) );
            }

            /**
            * @brief CRC-8/WCDMA with poly = 0x9B init = 0x00 refIn = true 
            * refOut = true  xorOut= 0x00 
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint8_t crc8_WCDMA( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint8_t>( generic( CRC8, pData, length, 0x9BU, 0U, true, true ) );
            }

            /**
            * @brief CRC-16/CCITT-FALSE with poly = 0x1021 init = 0xFFFF 
            * refIn = false refOut = false  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_CCITT_FALSE( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xFFFFU ) );
            }

            /**
            * @brief CRC-16/ARC with poly = 0x8005 init = 0x0000 refIn = true 
            * refOut = true  xorOut= 0x0000 
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_ARC( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0U, true, true ) );
            }

            /**
            * @brief CRC-16/AUG-CCITT with poly = 0x1021 init = 0x1D0F refIn = false
            * refOut = false  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_AUG_CCITT( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0x1D0FU ) );
            }

            /**
            * @brief CRC-16/BUYPASS with poly = 0x8005 init = 0x0000 refIn = false
            * refOut = false  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_BUYPASS( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U ) );
            }

            /**
            * @brief CRC-16/CDMA2000 with poly = 0xC867 init = 0xFFFF refIn = false
            * refOut = false  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_CDMA2000( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0xC867U, 0xFFFFU ) );
            }

            /**
            * @brief CRC-16/DDS-110 with poly = 0x8005 init = 0x800D refIn = false
            * refOut = false  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_DDS_110( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0x800DU ) );
            }

            /**
            * @brief CRC-16/DECT-R with poly = 0x0589 init = 0x0000 refIn = false
            * refOut = false  xorOut= 0x0001 
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_DECT_R( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x0589U, 0U, false, false, 0x0001U ) );
            }

            /**
            * @brief CRC-16/DECT-X with poly = 0x0589 init = 0x0000 refIn = false
            * refOut = false  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_DECT_X( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x0589U ) );
            }

            /**
            * @brief CRC-16/DNP with poly = 0x3D98 init = 0x0000 refIn = true 
            * refOut = true xorOut= 0xFFFF
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_DNP( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x3D65U, 0U, true, true, 0xFFFFU ) );
            }

            /**
            * @brief CRC-16/EN-13757 with poly = 0x3D65 init = 0x0000 refIn = false
            * refOut = false  xorOut= 0xFFFF
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_EN_13757( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x3D65U, 0U, false, false, 0xFFFFU ) );
            }

            /**
            * @brief CRC-16/GENIBUS with poly = 0x1021 init = 0xFFFF refIn = false 
            * refOut = false  xorOut= 0xFFFF
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_GENIBUS( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xFFFFU, false, false, 0xFFFFU ) );
            }

            /**
            * @brief CRC-16/MAXIM with poly = 0x8005 init = 0x0000 refIn = true 
            * refOut = true xorOut= 0xFFFF
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_MAXIM( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0U, true, true, 0xFFFFU ) );
            }

            /**
            * @brief CRC-16/MCRF4XX with poly = 0x1021 init = 0xFFFF refIn = true 
            * refOut = true xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_MCRF4XX( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xFFFFU, true, true ) );
            }

            /**
            * @brief CRC-16/RIELLO with poly = 0x1021 init = 0xB2AA refIn = true 
            * refOut = true  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_RIELLO( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xB2AAU, true, true ) );
            }

            /**
            * @brief CRC-16/DECT-X with poly = 0x8bb7 init = 0x0000 refIn = false 
            * refOut = false  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_T10_DIF( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8BB7U ) );
            }

            /**
            * @brief CRC-16/TELEDISK with poly = 0xA097 init = 0x0000 refIn = false 
            * refOut = false  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_TELEDISK( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0xA097U ) );
            }

            /**
            * @brief CRC-16/TMS37157 with poly = 0x1021 init = 0x89EC refIn = true 
            * refOut = true  xorOut= 0x000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_TMS37157( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0x89ECU, true, true ) );
            }

            /**
            * @brief CRC-16/USB with poly = 0x8005 init = 0xFFFF refIn = true 
            * refOut = true  xorOut= 0xFFFF
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_USB( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0xFFFFU, true, true, 0xFFFFU ) );
            }

            /**
            * @brief CRC-A with poly = 0x1021 init = 0xC6C6 refIn = true 
            * refOut = true  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_A( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xc6c6U, true, true ) );
            }

            /**
            * @brief CRC-16/KERMIT with poly = 0x1021 init = 0x0000 refIn = true 
            * refOut = true  xorOut= 0x000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_KERMIT( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0U, true, true ) );
            }

            /**
            * @brief CRC-16/MODBUS with poly = 0x8005 init = 0xFFFF refIn = true 
            * refOut = true  xorOut= 0x0000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_MODBUS( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x8005U, 0xFFFFU, true, true ) );
            }

            /**
            * @brief CRC-16/X-25 with poly = 0x1021 init = 0xFFFF refIn = true 
            * refOut = true  xorOut= 0xFFFF
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_X_25( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U, 0xFFFFU, true, true, 0xFFFFU ) );
            }

            /**
            * @brief CRC-16/XMODEM with poly = 0x1021 init = 0x0000 refIn = false 
            * refOut = false  xorOut= 0x0000 
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint16_t crc16_XMODEM( const void * const pData, const size_t length ) noexcept
            {
                return static_cast<uint16_t>( generic( CRC16, pData, length, 0x1021U ) );
            }

            /**
            * @brief CRC-32 with poly = 0x04C11DB7 init = 0xFFFFFFFF refIn = true 
            * refOut = true xorOut= 0xFFFFFFFF
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint32_t crc32( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0xFFFFFFFFU, true, true, 0xFFFFFFFFU  );
            }

            /**
            * @brief CRC-32/BZIP2 with poly = 0x04C11DB7 init = 0xFFFFFFFF 
            * refIn = false refOut = false xorOut= 0xFFFFFFFF 
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint32_t crc32_BZIP2( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0xFFFFFFFFU, false, false, 0xFFFFFFFFU  );
            }

            /**
            * @brief CRC-32C with poly = 0x1EDC6F41 init = 0xFFFFFFFF refIn = true 
            * refOut = true  xorOut= 0xFFFFFFFF
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint32_t crc32_C( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x1EDC6F41U, 0xFFFFFFFFU, true, true, 0xFFFFFFFFU  );
            }

            /**
            * @brief CRC-32D with poly = 0xA833982B init = 0xFFFFFFFF refIn = true 
            * refOut = true  xorOut= 0xFFFFFFFF 
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint32_t crc32_D( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0xA833982BU, 0xFFFFFFFFU, true, true, 0xFFFFFFFFU  );
            }

            /**
            * @brief CRC-32/JAMCRC with poly = 0x04C11DB7 init = 0xFFFFFFFF 
            * refIn = true refOut = true  xorOut= 0x00000000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint32_t crc32_JAMCRC( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0xFFFFFFFFU, true, true );
            }

            /**
            * @brief CRC-32/MPEG2 with poly = 0x04C11DB7 init = 0xFFFFFFFF 
            * refIn = false refOut = false  xorOut= 0x00000000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint32_t crc32_MPEG2( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0xFFFFFFFFU );
            }

            /**
            * @brief CRC-32/POSIX with poly = 0x04C11DB7 init = 0x00000000 
            * refIn = false refOut = false  xorOut= 0xFFFFFFFF 
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint32_t crc32_POSIX( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x04C11DB7U, 0U, false, false, 0xFFFFFFFFU  );
            }

            /**
            * @brief CRC-32Q with poly = 0x814141AB init = 0x00000000 
            * refIn = false refOut = false  xorOut= 0x00000000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint32_t crc32_Q( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x814141ABU );
            }

            /**
            * @brief CRC-32/XFER with poly = 0x000000AF init = 0x00000000 
            * refIn = false refOut = false  xorOut= 0x00000000
            * @param[in] pData A pointer to the block of data.
            * @param[in] length The number of bytes in @a data.
            * @return The CRC value for @a data.
             */
            static inline uint32_t crc32_XFER( const void * const pData, const size_t length ) noexcept
            {
                return generic( CRC32, pData, length, 0x000000AFU );
            }
    };

    /** @}*/

}

#endif /*QLIBS_CRC*/