/*
 * ViabilityTest.cpp
 *
 *  Created on: Thursday, September 21st, 2017
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

#include <xassert/XAssert.h>
#include <xutils/Log.h>
#include <xutils/Utils.h>
#include <xutils/Timer.h>

#include <bilinear/Library.h>
#include <bilinear/FastMultExp.h>
#include <bilinear/AppMain.h>

using namespace std;
using namespace Bilinear;

template<class GT>
void testAddMult();

template<class GT>
void testFastMultExp();

void testBNTbits();

void testFastModulo(const BNT& fieldOrder);

void testInvertModPrimeHelper(const BNT& a, const BNT& fieldOrder);
void testInvertModP(const BNT& fieldOrder);

int BilinearAppMain(const Library& lib, const std::vector<std::string>& args) {
    (void)lib;
    (void)args;

    loginfo << "RELIC Type 1 paring: " << pc_map_is_type1() << endl;
    loginfo << "RELIC Type 3 paring: " << pc_map_is_type3() << endl;

    unsigned int seed = static_cast<unsigned int>(time(NULL));
    loginfo << "Randomness seed passed to srand(): " << seed << endl;
    srand(seed);
    loginfo << endl;

    BNT fieldOrder = lib.getGroupOrder();

    int numIters = 100;
    logdbg << "Testing fast modular reduction..." << endl;
    for(int i = 0; i < numIters; i++) {
        testFastModulo(fieldOrder);
    }

    logdbg << "Testing accumulated multiplication and addition in G1..." << endl;
    for(int i = 0; i < numIters; i++) {
        testAddMult<G1T>();
    }

    logdbg << "Testing accumulated multiplication and addition in G2..." << endl;
    for(int i = 0; i < numIters; i++) {
        testAddMult<G2T>();
    }

    logdbg << "Testing RELIC's and BNT's getBit() ..." << endl;
    for(int i = 0; i < numIters; i++) {
        testBNTbits();
    }

    logdbg << "Testing RELIC's BNT::invertModPrime() ..." << endl;
    testInvertModPrimeHelper(BNT::One(), fieldOrder);
    testInvertModPrimeHelper(fieldOrder - BNT(1), fieldOrder);
    testInvertModPrimeHelper(fieldOrder - BNT(2), fieldOrder);
    for(int i = 0; i < numIters; i++) {
        testInvertModP(fieldOrder);
    }

    // NOTE: Now we are pretending G1 and G2 are groups using multiplicative notation
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

void testBNTbits() {
    BNT n(rand() % 544232);
    //BNT n(4);
    int a = 0;
    int bits = n.getBits();
    testAssertLessThanOrEqual(static_cast<size_t>(bits), (sizeof a) * 8);

    //for(int i = bits - 1; i >= 0; i--) {
    for(int i = 0; i < bits; i++) {
        bool bit = n.getBit(i);
        logtrace << "Bit " << i << " of " << n << " = " << static_cast<int>(bit) << endl;

        if(bit) {
            a += (1 << i);
        }
    }

    logtrace << n << " =?= " << a << endl;

    testAssertEqual(static_cast<int>(n.toDigit()), a);
}

template<class GT>
void testAddMult() {
    GT x, y;
    BNT e1, e2;
    BNT order;
    // Recall: G1, G2 and GT all have the same order
    g1_get_ord(order);
    x.Random();
    y.Random();
    e1.RandomMod(order);
    e2.RandomMod(order);

    GT r1, r2, r3;

    r1 = GT::Add(x, y);
    r2 = x; r2.Add(y);
    testAssertEqual(r1, r2);

    r1 = GT::Times(x, e1);
    r2 = x; r2.Times(e1);
    testAssertEqual(r1, r2);

    r1 = GT::Double(x);
    r2 = x; r2.Double();
    r3 = GT::Add(x, x);
    testAssertEqual(r1, r2);
    testAssertEqual(r1, r3);

    r1 = GT::TimesTwice(x, e1, y, e2);
    r2 = GT::Add(GT::Times(x, e1), GT::Times(y, e2));
    testAssertEqual(r1, r2);
}

template<class GT>
void testFastMultExp() {
    GT r1, r2;
    int n = 10 + (rand() % 2);
    int k = 5;
    // For fast multiple exponentiation, we need to know the max number of bits in an exponent
    int maxBits = Library::Get().getGroupOrderNumBits();
    //int maxBits = 256;
    //logdbg << "Max bits: " << maxBits << endl;

    std::vector<size_t> s;
    Utils::randomSubset(s, n, k);

    std::vector<GT> a;
    std::vector<BNT> e;
    a.resize(static_cast<size_t>(n));
    e.resize(static_cast<size_t>(n));

    for(size_t i : s) {
        a[i].Random();
        e[i].RandomMod(Library::Get().getGroupOrder());
    }

    testAssertEqual(r1, GT::Identity());
    testAssertEqual(r2, GT::Identity());

    // Fast way
    r1 = fastMultExp<GT>(s, a, e, maxBits);

    // Slow way
    for(size_t i : s) {
        GT& base = a[i];
        BNT& exp = e[i];

        GT pow = GT::Times(base, exp);
        r2.Add(pow);
    }

    // Same way?
    testAssertEqual(r1, r2);
}

void testFastModulo(const BNT& fieldOrder) {
    BNT u = BNT::FastModuloPre(fieldOrder);
    BNT a;
    // Pick random a > fieldOrder so we can reduce it
    a.RandomMod(fieldOrder);
    a.Times(a + BNT(3) + fieldOrder);

    BNT b = a;
    BNT pmers = a, monty = a, barrt = a;

    a.SlowModulo(fieldOrder);
    b.FastModulo(fieldOrder, u);
    testAssertEqual(a, b);

    u = BNT::FastModuloPrePmers(fieldOrder);
    pmers.FastModuloPmers(fieldOrder, u);
    testAssertEqual(a, pmers);

    u = BNT::FastModuloPreMonty(fieldOrder);
    monty.FastModuloMonty(fieldOrder, u);
    testAssertEqual(a, monty);

    u = BNT::FastModuloPreBarrett(fieldOrder);
    barrt.FastModuloBarrett(fieldOrder, u);
    testAssertEqual(a, barrt);
}

void testInvertModPrimeHelper(const BNT& a, const BNT& fieldOrder) {
    BNT b = a.invertModPrime(fieldOrder); 

    testAssertEqual((a*b).SlowModulo(fieldOrder), BNT::One());
}

void testInvertModP(const BNT& fieldOrder) {
    BNT a;
    a.RandomMod(fieldOrder);
    testInvertModPrimeHelper(a, fieldOrder);
}
