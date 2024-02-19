#include "qlibs.h"

using namespace qlibs;

struct thing{
    int a;
    float b;
};


bool operator<(const thing& lhs, const thing& rhs) {
    return ( lhs.a < rhs. a);
}


int main() {
    thing things[] = {
        {1,0},
        {-2,5},
        {8,4},
        {-2, 2},
    };
    algorithm::sort( things );
    algorithm::reverse( things );
    algorithm::rotate( things );
/*
    mat<2,2> x( 
                1.0f, 2.0f,
                3.0f, 4.0f
              );
    mat<2,1> y(
                1.0f,
                3.0f
              );
    mat<2,1> z( 
                2.0f,
                2.0f
              );

    mat<2,1> result = 4.5 - z - 2.0f*x*y - 3;
    x*=x;
    result.display();
    x.display();
  */
  return 0;
}
