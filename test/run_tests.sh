#! /bin/bash


if ! [ -x "$(command -v pytest)" ]; then
    echo 'Error: 'pytest' is not installed.' >&2
    exit 1
fi


# Get path to PFASSE root directory
PFASSTDIR="$(dirname $0)/.."

# Convert to absolute path
PFASSTDIR="$(realpath $PFASSTDIR)"

TESTDIR="$PFASSTDIR/test"

echo "Using PFASSTDIR=$PFASSTDIR"
echo "Using TESTDIR=$TESTDIR"


# Change to test directory
cd "$TESTDIR"
test -h "pf" || ln -s "$PFASSTDIR/pf"

cd "$TESTDIR/magpicard"
test -h "pf" || ln -s "$PFASSTDIR/pf"

cd "$TESTDIR/imk"
test -h "pf" || ln -s "$PFASSTDIR/pf"



cd "$PFASSTDIR"
OMP=y make -j || exit 1

cd "$TESTDIR/magpicard"
OMP=y make -j || exit 1

cd "$TESTDIR/imk"
OMP=y make -j || exit 1

cd "$TESTDIR/nagumo"
make -j || exit 1

cd "$TESTDIR/adv_diff_fft/1d"
make -j || exit 1

cd "$TESTDIR/adv_diff_fft/2d"
make -j


cd "$PFASSTDIR"
pytest
