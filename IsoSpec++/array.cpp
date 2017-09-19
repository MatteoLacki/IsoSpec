#include <cstdlib>
#include "array.h"

struct DoubleArray
{
    double *array;
    int used;
    int size;
};

void DoubleArrayConstructor(DoubleArray *a, int initialSize)
{
    a->array = (double *) malloc(initialSize * sizeof(double));
    a->used = 0;
    a->size = initialSize;
}

void DoubleArrayDestructor(DoubleArray *a, bool free_array = true)
{
    if( free_array ){ free(a->array); }
    a->array = NULL;
    a->used = a->size = 0;
}

void DoubleArrayPush(DoubleArray *a, double element)
{
    if (a->used == a->size)
    {
        a->size *= 2;
        a->array = (double *) realloc(a->array, a->size * sizeof(double));
    }
    a->array[a->used] = element;
    a->used++;
}
