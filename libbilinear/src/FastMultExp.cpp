/*
 * FastMultExp.cpp
 *
 *  Created on: March 24th, 2018
 *      Author: Alin Tomescu <alinush@mit.edu>
 */

#include <vector>

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
    const std::vector<GT>& a, 
    const std::vector<BNT>& e, int maxBits)
{
  (void) a;
  (void) e;
  (void) maxBits;
  return GT::Identity();
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

