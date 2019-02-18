/*
 * FastMultExp.h
 *
 *  Created on: March 24th, 2018
 *      Author: Alin Tomescu <alinush@mit.edu>
 */

#pragma once

#include "bilinear/Groups.h"

#include <vector>

using Bilinear::BNT;

/**
 * Computes r = \prod_{i \in S} { a[i]^e[i] } faster than the naive method using the result from
 * "Fast Batch Verification for Modular Exponentiation and Digital Signatures" by Bellare et al.
 */
template<class GT>
GT fastMultExp(
    const std::vector<size_t>& S, const std::vector<GT>& a, 
    const std::vector<BNT>& e, int maxBits);

template<class GT>
GT fastMultExp(
    const std::vector<GT>& a, 
    const std::vector<BNT>& e, int maxBits);

template<class GT>
GT fastMultExp_de_Rooij(
    const std::vector<GT>& a, 
    const std::vector<BNT>& e, int maxBits);
