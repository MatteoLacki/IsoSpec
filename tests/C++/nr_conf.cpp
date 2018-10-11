#include <iostream>
#include <string>
#include <cassert>

#include "../../IsoSpec++/unity-build.cpp"

int main()
{
  std::string formula = "C100";

  int tabSize = 1000;
  int hashSize = 1000;
  double threshold = 0.01;
  bool absolute = false;

  threshold = 1e-2;
  {
    Iso* iso = new Iso(formula.c_str());
    IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 
    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, true, true, false, true); 
    int64_t size = tabulator->confs_no();
    std::cout << size << std::endl;

    delete iso;
    delete generator;
    delete tabulator;
  }


  threshold = 1e-200;
  {
    Iso* iso = new Iso(formula.c_str());
    IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 
    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, true, true, false, true); 
    int64_t size = tabulator->confs_no();
    std::cout << size << std::endl;

    delete iso;
    delete generator;
    delete tabulator;
  }

  formula = "C520H817N139O147S8";

  threshold = 1e-10;
  {
    Iso* iso = new Iso(formula.c_str());
    IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 
    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, true, true, false, true); 
    int64_t size = tabulator->confs_no();
    std::cout << size << std::endl;

    delete iso;
    delete generator;
    delete tabulator;
  }

  threshold = 1e-50;
  {
    Iso* iso = new Iso(formula.c_str());
    IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 
    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, true, true, false, true); 
    int64_t size = tabulator->confs_no();
    std::cout << size << std::endl;

    delete iso;
    delete generator;
    delete tabulator;
  }

  threshold = 1e-100;
  {
    Iso* iso = new Iso(formula.c_str());
    IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 
    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, true, true, false, true); 
    int64_t size = tabulator->confs_no();
    std::cout << size << std::endl;

    delete iso;
    delete generator;
    delete tabulator;
  }

#if 0
  threshold = 1e-200;
  {
    Iso* iso = new Iso(formula.c_str());
    IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 
    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, true, true, false, true); 
    int64_t size = tabulator->confs_no();
    std::cout << size << std::endl;

    delete iso;
    delete generator;
    delete tabulator;
  }

  threshold = 1e-300;
  {
    Iso* iso = new Iso(formula.c_str());
    IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 
    Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, true, true, false, true); 
    int64_t size = tabulator->confs_no();
    std::cout << size << std::endl;

    delete iso;
    delete generator;
    delete tabulator;
  }
#endif
}
