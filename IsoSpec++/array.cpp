extern "C"{

#include <stdio.h>
#include <stdlib.h>
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

void DoubleArrayDestructor(DoubleArray *a)
{
    free(a->array);
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


}
