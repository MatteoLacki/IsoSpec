// Regression test for: fasta.h uses memset() without including <cstring>.
//
// This TU includes ONLY fasta.h.  It must compile cleanly.  Today it does
// not: parse_fasta's inline definition references memset, and that symbol
// is not declared by any header fasta.h pulls in.
//
// Build:
//   clang++ -std=c++17 -I../../src/IsoSpec++ -c fasta_header_isolation.cpp
//
// Expected behaviour after fix: compiles to an object file without error.
// Expected behaviour today: fails with "use of undeclared identifier 'memset'".

#include "fasta.h"

using namespace IsoSpec;

// Reference the inline so the TU actually instantiates it.
int probe(const char* s, int out[6])
{
    parse_fasta(s, out);
    return out[0];
}
