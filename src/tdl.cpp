#include <include/tdl.hpp>

using namespace qlibs;

/*============================================================================*/
void tdl::setup( real_t * const area,
                 const size_t n,
                 const real_t initVal ) noexcept
{
    itemCount = n;
    head = area;
    flush( initVal );
}
/*============================================================================*/
void tdl::flush( const real_t initVal ) noexcept
{
    tail = &head[ itemCount ];
    wr = head;
    /*cstat -CERT-INT30-C_a*/
    rd = &head[ itemCount - 1U ];
    /*cstat +CERT-INT30-C_a*/
    for ( size_t i = 0U ; i < itemCount ; ++i ) { /*initialize the queue*/
        insertNewest( initVal ); /*the queue should be always full*/
    }
}
/*============================================================================*/
void tdl::insertNewest( const real_t sample ) noexcept
{
    wr[ 0 ] = sample;
    wr++;
    if ( wr >= tail ) {
        wr = head;
    }
}
/*============================================================================*/
void tdl::removeOldest( void ) noexcept
{
    if ( ++rd >= tail ) {
        rd = head;
    }
}
/*============================================================================*/
real_t tdl::getOldest( void ) const noexcept
{
    return ( ( rd + 1U ) >= tail ) ? head[ 0 ] : rd[ 1 ];
}
/*============================================================================*/
real_t tdl::getRecent( void ) const noexcept
{
    return rd[ 0 ];
}
/*============================================================================*/
real_t tdl::getAtIndex( const size_t i ) const noexcept
{
    /*cstat -MISRAC++2008-5-0-3*/
    return ( ( wr >= rd ) && ( ( head + i ) >= wr ) ) ? rd[ itemCount - i ]
                                                      : *( rd - i );
    /*cstat +MISRAC++2008-5-0-3*/
}
/*============================================================================*/
void tdl::insertSample( const real_t sample ) noexcept
{
    removeOldest();
    insertNewest( sample );
}
/*============================================================================*/
const real_t& tdl::operator[]( int index ) noexcept
{
    const size_t i = static_cast<size_t>( index );
    /*cstat -MISRAC++2008-6-6-5*/
    if ( ( index >= 0 ) && ( i < itemCount ) ) {
        /*cstat -MISRAC++2008-5-0-3*/
        return ( ( wr >= rd ) && ( ( head + i ) >= wr ) ) ? rd[ itemCount - i ]
                                                          : *( rd - i );
        /*cstat +MISRAC++2008-5-0-3*/
    }
    else if ( -1 == index ) {
        return ( ( rd + 1U ) >= tail ) ? head[ 0 ] : rd[ 1 ];
    }
    else {
        return undefined;
    }
    /*cstat +MISRAC++2008-6-6-5*/
}
/*============================================================================*/
