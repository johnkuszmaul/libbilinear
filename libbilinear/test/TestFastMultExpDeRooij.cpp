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

int BilinearAppMain(const Library& lib, const std::vector<std::string>& args) {
    (void)args;
    (void)lib;

    // TODO: test the multiexp
    fastMultExp_de_Rooij<G1T>(std::vector<G1T>(), std::vector<BNT>(), 0);

    return 0;
}


