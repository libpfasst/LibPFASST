#! /bin/bash


if ! [ -x "$(command -v pytest)" ]; then
    echo 'Error: 'pytest' is not installed.' >&2
    exit 1
fi


# Get path to PFASS root directory
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

echo "Building libpfasst"
cd "$PFASSTDIR"
make clean
make -j || exit 1

echo "Building magpicard"
cd "$TESTDIR/magpicard"
make clean
make -j || exit 1

echo "Building imk"
cd "$TESTDIR/imk"
make clean
OMP=y make -j || exit 1

echo "Building nagumo"
cd "$TESTDIR/nagumo"
make clean
make -j || exit 1

echo "Building EXP 1d"
cd "$TESTDIR/EXP_adv_diff_fft/1d"
make clean
make -j || exit 1

echo "Building ad 1d"
cd "$TESTDIR/adv_diff_fft/1d"
make clean
make -j || exit 1

echo "Building ad 2d"
cd "$TESTDIR/adv_diff_fft/2d"
make clean
make -j


cd "$PFASSTDIR"
pytest
