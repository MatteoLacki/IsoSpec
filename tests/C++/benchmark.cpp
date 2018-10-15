#include <iostream>
#include <string>
#include <cassert>

#include <ctime>

#include "../../IsoSpec++/unity-build.cpp"

int main()
{
  std::string formula = "C100";
  formula = "C520H817N139O147S8";

  int tabSize = 1000;
  int hashSize = 1000;
  double threshold = 0.01;
  bool absolute = false;

  formula = "C520H817N139O147S8";
  threshold = 1e-95; // 2,000M entries -> ca 31 GB of RAM for both arrays
  threshold = 1e-75; // 755M entries -> ca. 11.8 GB of RAM for both arrays
  threshold = 1e-60; // 295M entries -> ca 4.6 GB of RAM for both arrays
  threshold = 1e-50; // 135M entries -> ca 2.1 GB of RAM for both arrays
  threshold = 1e-45; // 85M entries -> ca 1.3 GB of RAM for both arrays

  threshold = 1e-95; 

  int64_t tsize;

  #if (false)
  {
    clock_t begin = std::clock();
    Iso iso(formula.c_str());
    IsoThresholdGeneratorCntr generator (std::move(iso), threshold, absolute, tabSize, hashSize); 
    int64_t size = 0;
    while ((generator.advanceToNextConfiguration())) size++;
    std::cout << size << std::endl;
    clock_t end = std::clock();
    std::cout << "Using IsoThresholdGeneratorCntr, it took " <<   double(end - begin) / CLOCKS_PER_SEC << "s "<< std::endl;
    tsize = size;
  }
  #endif
  // fast
  #if (false)
  {
    clock_t begin = std::clock();
    Iso iso(formula.c_str());
    IsoThresholdGeneratorFast generator (std::move(iso), threshold, absolute, tabSize, hashSize); 
    int64_t size = 0;
    while (generator.advanceToNextConfiguration()) size++; // std::cout << size << std::endl;}
    std::cout << size << std::endl;
    clock_t end = std::clock();
    std::cout << "Using IsoThresholdGeneratorFast, it took " <<   double(end - begin) / CLOCKS_PER_SEC << "s "<< std::endl;
    tsize = size;
  }
  #endif

  // legacy
  if (true)
  {
    clock_t begin = std::clock();
    Iso iso(formula.c_str());
    IsoThresholdGenerator generator (std::move(iso), threshold, absolute, tabSize, hashSize); 
    int64_t size = 0;
    while (generator.advanceToNextConfiguration()) size++;
    std::cout << size << std::endl;
    clock_t end = std::clock();
    std::cout << "Using new, it took " <<   double(end - begin) / CLOCKS_PER_SEC << "s "<< std::endl;
    tsize = size;
  }

  // fast, calculating a (meaningless) sum of masses
  #if (false)
  {
    clock_t begin = std::clock();
    Iso iso(formula.c_str());
    IsoThresholdGeneratorFast generator (std::move(iso), threshold, absolute, tabSize, hashSize);
    double total = 0.0;
    while (generator.advanceToNextConfiguration()) total += generator.mass();
    std::cout << total << std::endl;
    clock_t end = std::clock();
    std::cout << "Sum of masses: using IsoThresholdGeneratorFast, it took " <<   double(end - begin) / CLOCKS_PER_SEC << "s "<< std::endl;
  }
  #endif

  // legacy, 
  #if (true)
  {
    clock_t begin = std::clock();
    Iso iso(formula.c_str());
    IsoThresholdGenerator generator (std::move(iso), threshold, absolute, tabSize, hashSize); 
    double total = 0.0;
    while (generator.advanceToNextConfiguration()) total += generator.lprob();
    std::cout << total << std::endl;
    clock_t end = std::clock();
    std::cout << "sum of lprobs: Using new, it took " <<   double(end - begin) / CLOCKS_PER_SEC << "s "<< std::endl;
  }
  #endif

  // RAM-based
  if (tsize <= 3000000000)
  {
    double* T = new double[tsize];
    // Fill it with junk
    for(int ii=0; ii<tsize; ii++)
        T[ii] = 17.0456;
    double sum = 0.0;
    clock_t begin = std::clock();
    for(int ii=0; ii<tsize; ii++)
        sum += T[ii];
    clock_t end = std::clock();
    delete[] T;
    std::cout << sum << std::endl;
    std::cout << "Memory sum: it took " <<   double(end - begin) / CLOCKS_PER_SEC << "s "<< std::endl;
  }


  // verify that the results are the same
  #if(fasle)
  {
    IsoThresholdGenerator generator_old(formula.c_str(), threshold, absolute, tabSize, hashSize);
    IsoThresholdGeneratorFast generator_new(formula.c_str(), threshold, absolute, tabSize, hashSize);
    int confspace_size = generator_old.getAllDim();
    int confspace_old[confspace_size];
    int confspace_new[confspace_size];

    while(generator_new.advanceToNextConfiguration())
    {
        assert(generator_old.advanceToNextConfiguration());
        assert(generator_new.eprob() == generator_old.eprob());
        assert(generator_new.lprob() == generator_old.lprob());
        assert(generator_new.mass()  == generator_old.mass());
        generator_old.get_conf_signature(confspace_old);
        generator_new.get_conf_signature(confspace_new);
        assert(0 == memcmp(confspace_old, confspace_new, confspace_size*sizeof(int)));
    }
    assert(!generator_old.advanceToNextConfiguration());
    std::cout << "Results are the same for old and new!" << std::endl;

  }
  #endif
}
