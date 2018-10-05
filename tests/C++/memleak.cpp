#include <iostream>
#include <string>
#include <cassert>

#include "isoSpec++.h"
#include "tabulator.cpp"

int main()
{
  std::string formula = "C520H817N139O147";

  int tabSize = 1000;
  int hashSize = 1000;
  double threshold = 0.01;
  bool absolute = false;

  // Do some stress testing of the library...
  int sum = 0;
  for (size_t k = 0; k < 5e6; k++)
  {

    // Run threshold with 1e-2
    {
      Iso* iso = new Iso(formula.c_str());
      IsoThresholdGenerator* generator = new IsoThresholdGenerator(std::move(*iso), threshold, absolute, tabSize, hashSize); 
      Tabulator<IsoThresholdGenerator>* tabulator = new Tabulator<IsoThresholdGenerator>(generator, true, true, false, true); 
      int size = tabulator->confs_no();
      sum += size;
      // std::cout << size << std::endl;

      delete generator;
      delete tabulator;
    }


  }
  // with debug on:
  // 1e4 :   1.18user 0.10system 0:01.28elapsed 99%CPU (0avgtext+0avgdata 287836maxresident)k
  // 1e5 :  11.70user 0.14system 0:11.95elapsed 99%CPU (0avgtext+0avgdata 295528maxresident)k
  // 5e5 :  58.12user 0.33system 0:58.57elapsed 99%CPU (0avgtext+0avgdata 337772maxresident)k
  // 1e6 : 117.65user 0.20system 1:57.86elapsed 99%CPU (0avgtext+0avgdata 392104maxresident)k

  // gcc-opt
  // 1e5 :   1.66user 0.00system 0:01.66elapsed 99%CPU (0avgtext+0avgdata  9884maxresident)k
  // 1e6 :  16.21user 0.00system 0:16.22elapsed 99%CPU (0avgtext+0avgdata  66120maxresident)k
  // 5e6 :  86.24user 0.17system 1:26.41elapsed 99%CPU (0avgtext+0avgdata 316224maxresident)k
  //
  // clang-opt
  // 1e5 :   1.66user 0.00system 0:01.66elapsed 99%CPU (0avgtext+0avgdata  10060maxresident)k
  // 1e6 :  16.79user 0.02system 0:16.81elapsed 99%CPU (0avgtext+0avgdata  66192maxresident)k
  // 5e6 : 103.04user 0.17system 1:43.24elapsed 99%CPU (0avgtext+0avgdata 316348maxresident)k
  //
  //

}


