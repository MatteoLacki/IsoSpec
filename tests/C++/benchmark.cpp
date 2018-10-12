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

  clock_t begin = std::clock();
  formula = "C520H817N139O147S8";
  threshold = 1e-95; // 2,000M entries -> ca 31 GB of RAM for both arrays
  threshold = 1e-75; // 755M entries -> ca. 11.8 GB of RAM for both arrays
  threshold = 1e-60; // 295M entries -> ca 4.6 GB of RAM for both arrays
  threshold = 1e-50; // 135M entries -> ca 2.1 GB of RAM for both arrays
  threshold = 1e-45; // 85M entries -> ca 1.3 GB of RAM for both arrays

  threshold = 1e-95; 

  if (true)
  {
    clock_t begin = std::clock();
    Iso iso(formula.c_str());
    IsoThresholdGeneratorCntr generator (std::move(iso), threshold, absolute, tabSize, hashSize); 
    int64_t size = 0;
    while ((generator.advanceToNextConfiguration())) size++;
    std::cout << size << std::endl;
    clock_t end = std::clock();
    std::cout << "Using IsoThresholdGeneratorCntr, it took " <<   double(end - begin) / CLOCKS_PER_SEC << "s "<< std::endl;
  }

  // fast
  if (true)
  {
    clock_t begin = std::clock();
    Iso iso(formula.c_str());
    IsoThresholdGeneratorFast generator (std::move(iso), threshold, absolute, tabSize, hashSize); 
    int64_t size = 0;
    while (generator.advanceToNextConfiguration()) size++; // std::cout << size << std::endl;}
    std::cout << size << std::endl;
    clock_t end = std::clock();
    std::cout << "Using IsoThresholdGeneratorFast, it took " <<   double(end - begin) / CLOCKS_PER_SEC << "s "<< std::endl;
  }


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
    std::cout << "Using legacy, it took " <<   double(end - begin) / CLOCKS_PER_SEC << "s "<< std::endl;
  }

}
