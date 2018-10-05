#include <iostream>
#include <string>
#include <cassert>

#include "isoSpec++.h"
#include "tabulator.cpp"

int main()
{
  std::string formula = "C6H12O6";
  formula = "H2O1";

  int tabSize = 1000;
  int hashSize = 1000;
  double threshold = 0.01;
  bool absolute = false;

#ifdef PREVENT_SEGFAULT
#else
  threshold = 1e-5;
  // Run threshold with 1e-5
  {
    Iso* iso = new Iso(formula.c_str());
    IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 
    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, true, true, false, true); 
    int size = tabulator->confs_no();
    // assert(size == 14);

    delete generator;
    delete tabulator;
  }
#endif

  // Run Layered
  {
    Iso* iso = new Iso(formula.c_str());

    double delta = -10;
    IsoLayeredGenerator* generator = new IsoLayeredGenerator(std::move(*iso), delta, tabSize, hashSize);
    Tabulator<IsoLayeredGenerator>* tabulator = new Tabulator<IsoLayeredGenerator>(generator, true, true, true, true); 

    int size = tabulator->confs_no();

    // assert(size == 32);

    delete generator;
    delete tabulator;
  }

}


