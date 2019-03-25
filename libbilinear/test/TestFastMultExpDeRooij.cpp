/*
 * ElementSizesTest.cpp
 *
 *      Author: alinush
 */


#include "bilinear/Configuration.h"

#include <map>
#include <set>
#include <vector>
#include <string>
#include <cassert>
#include <memory>
#include <stdexcept>
#include <inttypes.h>

#include "xutils/Log.h"
#include "xutils/Utils.h"
#include "xutils/Timer.h"
#include "xassert/XAssert.h"

#include "bilinear/Library.h"
#include "bilinear/AppMain.h"
#include "bilinear/FastMultExp.h"

using namespace std;
using namespace Bilinear;

template<class GT>
void testFastMultExp();

int BilinearAppMain(const Library& lib, const std::vector<std::string>& args) {
  (void)args;
  (void)lib;

  // TODO: test the multiexp
  //    fastMultExp_de_Rooij<G1T>(std::vector<G1T>(), std::vector<BNT>(), 0);

  int numIters = 100;    

  logdbg << "Testing fast exponentiated multiplication in G1..." << endl;
  for(int i = 0; i < numIters; i++) {
    testFastMultExp<G1T>();
  }

  logdbg << "Testing fast exponentiated multiplication in G2..." << endl;
  for(int i = 0; i < numIters; i++) {
    testFastMultExp<G2T>();
  }

  return 0;
}

template<class GT>
void testFastMultExp() {
  GT r1, r2, r3;
  size_t n = 10 + static_cast<size_t>(rand() % 2);
  // For fast multiple exponentiation, we need to know the max number of bits in an exponent
  int maxBits = Library::Get().getGroupOrderNumBits();

  std::vector<GT> a;
  std::vector<BNT> e;
  a.resize(static_cast<size_t>(n) + 1);
  e.resize(static_cast<size_t>(n) + 1);

  for(size_t i = 0; i <= n; i++) {
    a[i].Random();
    e[i].RandomMod(Library::Get().getGroupOrder());
  }

  assertEqual(r1, GT::Identity());
  assertEqual(r2, GT::Identity());
  assertEqual(r3, GT::Identity());

  logdbg << "Testing Fast way" << std::endl;

  // Fast way
  r1 = fastMultExp<GT>(a, e, maxBits);

  logdbg << "Testing De Rooij way" << std::endl;    
    
  // De Rooij Way
  r3 = fastMultExp_de_Rooij<GT>(a, e, maxBits); // TODO: Not defined for (s, a, e, maxbits)

  logdbg << "Testing slow way" << std::endl;    

  // Slow way
  for(size_t i = 0; i <= n; i++) {
    GT& base = a[i];
    BNT& exp = e[i];

    GT pow = GT::Times(base, exp);
    r2.Add(pow);
  }

  // Same way?
  assertEqual(r1, r2);
  assertEqual(r1, r3);

  //  logdbg << "Fast way: " << r1 << std::endl;
  //  logdbg << "Slow way: " << r2 << std::endl;
  //  logdbg << "De Rooij way: " << r3 << std::endl;
}
