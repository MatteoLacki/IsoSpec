#ifndef __ARRAYC
#define __ARRAYC

typedef struct DoubleArray DoubleArray;

void DoubleArrayConstructor(DoubleArray *a, int initialSize);
void DoubleArrayDestructor(DoubleArray *a, bool free_array);
void DoubleArrayPush(DoubleArray *a, double element);

#endif /* __ARRAYC */
