/*
 * FastMultExp.cpp
 *
 *  Created on: March 24th, 2018
 *      Author: Alin Tomescu <alinush@mit.edu>
 */

#include <vector>
#include <algorithm>
#include <numeric>

#include "bilinear/FastMultExp.h"
#include "bilinear/Groups.h"

#include <xutils/Log.h>
#include <xutils/Timer.h>
#include <xassert/XAssert.h>

using std::endl;
using namespace Bilinear;

template<class GT>
GT fastMultExp(
    const std::vector<size_t>& S, const std::vector<GT>& a, 
    const std::vector<BNT>& e, int maxBits)
{
    GT r;
    assertEqual(r, GT::Identity());

    for(int j = maxBits - 1; j >= 0; j--) {
        r.Double();

        for(size_t idx : S) {
            assertLessThanOrEqual(e[idx].getBits(), maxBits);

            if(e[idx].getBit(j))
                r.Add(a[idx]);
        }
    }

    return r;
}

template<class GT>
GT fastMultExp(
    const std::vector<GT>& a, 
    const std::vector<BNT>& e, int maxBits)
{
    GT r;
    assertEqual(r, GT::Identity());

    for(int j = maxBits - 1; j >= 0; j--) {
        r.Double();

        for(size_t idx = 0; idx < a.size(); idx++) {
            assertLessThanOrEqual(e[idx].getBits(), maxBits);

            if(e[idx].getBit(j))
                r.Add(a[idx]);
        }
    }

    return r;
}

template<class GT>
GT fastMultExp_de_Rooij(
    const std::vector<GT>& bases, 
    const std::vector<BNT>& exps, int maxBits)
{
    CumulativeTimer all;
    CumulativeTimer div, exp, add;
    CumulativeTimer overhead;

    all.start();
    (void) maxBits;

    if(exps.size() != bases.size()) {
        throw std::runtime_error("the sizes of the bases and exponents must match");
    }
    if(exps.empty()) {
        throw std::runtime_error("give me some exponents and bases man");
    }

    overhead.start();
    std::vector<GT> b(bases.begin(), bases.end());
    std::vector<BNT> t(exps.begin(), exps.end());

    // build a max heap over the exponents
    std::vector<size_t> heap(t.size());
    auto comp = [&t](const size_t lhs, const size_t rhs) {
       return t[lhs] < t[rhs];
    };
    std::iota(heap.begin(), heap.end(), 0); // puts 0, 1, ..., e.size() - 1
    std::make_heap(heap.begin(), heap.end(), comp);
    overhead.end();

    //logdbg << "Built a max heap over the " << heap.size() << " exponent(s)" << std::endl;

    // get the index of the max exponent
    BNT q;
    BNT qAvg = BNT::Zero();
    size_t e_max;
    overhead.start();
    std::pop_heap(heap.begin(), heap.end(), comp);
    e_max = heap.back();
    heap.pop_back();
    overhead.end();
    size_t numLaps = 0;
    while(heap.size() > 0) {
        size_t e_next = heap[0]; // the index of the next max exponent

        if(t[e_next] == 0) {
            break;
        }

        assertGreaterThanOrEqual(t[e_max], t[e_next]);

        // First, q = t_max / t_next
        // Second, t_max = t_max mod t_next
        div.start();
        bn_div_rem(q, t[e_max], t[e_max], t[e_next]);
        div.end();

        // Third, b_next = (b_max ^ q) * b_next
        numLaps++;
        qAvg = qAvg + q;

        exp.start();
        GT interm;
        if(q == BNT::One()) {
            interm = b[e_max];
        } else {
            interm = GT::Times(b[e_max], q);
        }
        exp.end();
        add.start();
        b[e_next] = GT::Add(interm, b[e_next]);
        add.end();

        // push back the modified e_max
        if(t[e_max] != 0) {
            overhead.start();
            heap.push_back(e_max);
            std::push_heap(heap.begin(), heap.end(), comp);
            overhead.end();
        }
        
        // pop the next e_max
        overhead.start();
        std::pop_heap(heap.begin(), heap.end(), comp);
        e_max = heap.back();
        heap.pop_back();
        overhead.end();
    }
  
    GT res = GT::Times(b[e_max], t[e_max]);
    all.end();
    
    loginfo << "All: " << all.getTotal() << endl;
    loginfo << "OVR: " << overhead.getTotal() << endl;
    loginfo << "   fraction: " << static_cast<double>(overhead.getTotal())/static_cast<double>(all.getTotal()) * 100.0 << "%" << std::endl;
    loginfo << "Div: " << div.getTotal() << endl;
    loginfo << "Add: " << add.getTotal() << endl;
    loginfo << "Exp: " << exp.getTotal() << endl;
    loginfo << "numLaps: " << numLaps << ", avg Q: " << qAvg << endl;

    return res;
}

/**
 * Template instatiations, since we only use these with G1 and G2
 */
template G1T fastMultExp<G1T>(
    const std::vector<size_t>& S, const std::vector<G1T>& a,
    const std::vector<BNT>& e, int maxBits);
template G2T fastMultExp<G2T>(
    const std::vector<size_t>& S, const std::vector<G2T>& a, 
    const std::vector<BNT>& e, int maxBits);

template G1T fastMultExp<G1T>(
    const std::vector<G1T>& a,
    const std::vector<BNT>& e, int maxBits);
template G2T fastMultExp<G2T>(
    const std::vector<G2T>& a, 
    const std::vector<BNT>& e, int maxBits);

template G1T fastMultExp_de_Rooij<G1T>(
    const std::vector<G1T>& a,
    const std::vector<BNT>& e, int maxBits);
template G2T fastMultExp_de_Rooij<G2T>(
    const std::vector<G2T>& a, 
    const std::vector<BNT>& e, int maxBits);

