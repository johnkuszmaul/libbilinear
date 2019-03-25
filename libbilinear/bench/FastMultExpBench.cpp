/*
 * FastMultExpBench.cpp
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

#include <xutils/Log.h>
#include <xutils/Utils.h>
#include <xutils/Timer.h>
#include <xassert/XAssert.h>

using namespace std;

#include "bilinear/Library.h"
#include "bilinear/FastMultExp.h"
#include "bilinear/AppMain.h"

using namespace Bilinear;

template<class GT>
void benchFastMultExp(int numIters, int n, bool skipNaive = false);

#define FAST_MULT_COMPARE 1
#define FAST_MULT_LARGE 2

int BilinearAppMain(const Library& lib, const std::vector<std::string>& args) {
    (void)lib;
    (void)args;

    int type = FAST_MULT_COMPARE;
    int n = 512, numIters;

    if((args.size() > 1 && (args[1] == "-h" || args[1] == "--help")) || args.size() < 4) {
        cout << endl;
        cout << "Usage: " << args[0] << " <benchmark-type> <num-exps> <num-iters>" << endl;
        cout << endl;
        cout << "OPTIONS: " << endl;
        cout << "   <benchmark-type> can be:" << endl;
        cout << "     comparison  -- compares fast batch exponentiation method with naive one (default benchmark)" << endl;
        cout << "     large       -- just runs the fast batch exponentiation method" << endl;
        cout << "   <num-exps> is the number of exponentiations that will be benchmarked" << endl;
        cout << "   <num-iters> is the number of times a multi exponentiation will be measured" << endl;
        cout << endl;
        return 1;
    }

    if(args[1] == "comparison") {
        // already set default values above
    } else if (args[1] == "large") {
        type = FAST_MULT_LARGE;
        n = 1024*1024;
    } else {
        cout << "ERROR: Method must be either 'comparison' or 'large'" << endl;
        return 1;
    }

    n = std::stoi(args[2]);

    numIters = std::stoi(args[3]);

    //srand(seed);

    switch(type) {
        case FAST_MULT_COMPARE:
        {
            loginfo << "Benchmarking " << n << " exponentiations in G1 (fast batch + deRooij + naive)..." << endl;
            benchFastMultExp<G1T>(numIters, n);
            loginfo << endl;
            
            loginfo << "Benchmarking " << n << " exponentiations in G2 (fast batch + deRooij + naive)..." << endl;
            benchFastMultExp<G2T>(numIters, n);
            break;
        }
        case FAST_MULT_LARGE:
        {
            loginfo << "Benchmarking " << n << " exponentiations in G1 (fast batch only)..." << endl;
            benchFastMultExp<G1T>(numIters, n, true);
            loginfo << endl;

            loginfo << "Benchmarking " << n << " exponentiations in G2 (fast batch only)..." << endl;
            benchFastMultExp<G2T>(numIters, n, true);
            break;
        }
        default:
            throw std::logic_error("You did not correctly implement the default benchmarking path");
    }

    return 0;
}

template<class GT>
void benchFastMultExp(int numIters, int n, bool skipNaive) {
    GT r1, r2, r3;
    int maxBits = Library::Get().getGroupOrderNumBits();
    logdbg << "Max bits: " << maxBits << endl;

    // Store random group elements in a[] and their exponents in e[]
    std::vector<GT> a;
    std::vector<BNT> e;
    a.resize(static_cast<size_t>(n) + 1);
    e.resize(static_cast<size_t>(n) + 1);

    for(size_t i = 0; i < a.size(); i++) {
        a[i].Random();
        e[i].RandomMod(Library::Get().getGroupOrder());
    }

    assertEqual(r1, GT::Identity());
    assertEqual(r2, GT::Identity());
    assertEqual(r3, GT::Identity());

    // Slow way
    AveragingTimer t1("naiveExp    ");
    if(skipNaive == false) {
        for(int i = 0; i < numIters; i++) {
            t1.startLap();
            r2 = GT::Identity();
            for(size_t i = 0; i < a.size(); i++) {

                GT& base = a[i];
                BNT& exp = e[i];

                GT pow = GT::Times(base, exp);
                r2.Add(pow);
            }
            t1.endLap();
        }
    }

    // Fast way
    AveragingTimer t2("fastMultExp ");
    for(int i = 0; i < numIters; i++) {
        t2.startLap();
        r1 = fastMultExp<GT>(a, e, maxBits);
        t2.endLap();
    }

    // DeRooij way
    AveragingTimer t3("deRooijExp  ");
    for(int i = 0; i < numIters; i++) {
        t3.startLap();
        r3 = fastMultExp_de_Rooij<GT>(a, e, maxBits);
        t3.endLap();
    }

    if(skipNaive == false) {
        loginfo << t1 << endl;
    }
    loginfo << t2 << endl;
    loginfo << t3 << endl;

    if(skipNaive == false) {
        // Same way?
        if(r1 != r2) {
            throw std::runtime_error("Incorrect results returned by one of the implementations: naive or fast.");
        }   
        if(r1 != r3) {
            throw std::runtime_error("Incorrect results returned by one of the implementations: naive or deRooij.");
        }   
    }
}

