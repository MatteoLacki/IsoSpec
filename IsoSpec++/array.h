#ifndef __ARRAYC
#define __ARRAYC

extern "C"{

#include <stdlib.h>
typedef struct DoubleArray DoubleArray;

void DoubleArrayConstructor(DoubleArray *a, int initialSize);
void DoubleArrayDestructor(DoubleArray *a);
void DoubleArrayPush(DoubleArray *a, double element);

}
#endif /* __ARRAYC */
