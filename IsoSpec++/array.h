#ifndef __ARRAYC
#define __ARRAYC

extern "C"{

#include <stdlib.h>
typedef struct DoubleArray DoubleTable;

void ArrayConstructor(DoubleTable *a, size_t initialSize);
void ArrayDestructor(DoubleTable *a);
void ArrayPush(DoubleTable *a, double element);

}
#endif /* __ARRAYC */
