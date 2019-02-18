/*
 * FastMultExp.cpp
 *
 *  Created on: March 24th, 2018
 *      Author: Alin Tomescu <alinush@mit.edu>
 */

#include <vector>
#include <algorithm>

#include "bilinear/FastMultExp.h"
#include "bilinear/Groups.h"

#include "xutils/Log.h"
#include "xassert/XAssert.h"

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
    (void) maxBits;

    if(exps.size() != bases.size()) {
        throw std::runtime_error("the sizes of the bases and exponents must match");
    }
    if(exps.empty()) {
        throw std::runtime_error("give me some exponents and bases man");
    }

    std::vector<GT> b(bases.begin(), bases.end());
    std::vector<BNT> t(exps.begin(), exps.end());

    // build a max heap over the exponents
    std::vector<size_t> heap(t.size());
    auto comp = [&t](const size_t lhs, const size_t rhs) {
       return t[lhs] < t[rhs];
    };
    std::iota(heap.begin(), heap.end(), 0); // puts 0, 1, ..., e.size() - 1
    std::make_heap(heap.begin(), heap.end(), comp);

    // get the index of the max exponent
    size_t e_max;
    std::pop_heap(heap.begin(), heap.end(), comp);
    e_max = heap.back();
    heap.pop_back();
    while(heap.size() > 0) {
        size_t e_next = heap[0]; // the index of the next max exponent

        // TODO: optimize this into a single RELIC op if possible (check BNT division API in RELIC)
        // i.e., use https://github.com/relic-toolkit/relic/blob/master/include/relic_bn.h#L833
        t[e_max] = t[e_max] % t[e_next];
        BNT q = t[e_max];
        q.DivideBy(t[e_next]);

        b[e_next] = GT::Add(GT::Times(b[e_max], q), b[e_next]);

        // push back the modified e_max
        if(e_max != 0) {
            heap.push_back(e_max);
            std::push_heap(heap.begin(), heap.end(), comp);
        }

        // pop the next e_max
        std::pop_heap(heap.begin(), heap.end(), comp);
        e_max = heap.back();
        heap.pop_back();
    }
  
    return GT::Times(b[e_max], t[e_max]);
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

