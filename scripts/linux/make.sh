#!/bin/bash
set -e

scriptdir=$(cd $(dirname $0); pwd -P)
sourcedir=$(cd $scriptdir/../..; pwd -P)
. $scriptdir/shlibs/check-env.sh
. $scriptdir/shlibs/os.sh

builddir=$LIBBILINEAR_BUILD_DIR

echo "Making in '$builddir' ..."
echo
[ ! -d "$builddir" ] && $scriptdir/cmake.sh
cd "$builddir"
make $@ -j $NUM_CPUS
