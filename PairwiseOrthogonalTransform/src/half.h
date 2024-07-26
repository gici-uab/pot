#ifndef HALF_H
#define HALF_H

#include<stdint.h>
//#include<cstdint> when C++0x support is ready

uint32_t half_to_float( uint16_t h );
uint16_t half_from_float( uint32_t f );

#endif /* HALF_H */
