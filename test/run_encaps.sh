#! /bin/bash


#if ! [ -x "$(command -v pytest)" ]; then
#    echo 'Error: 'pytest' is not installed.' >&2
#    exit 1
#fi


# Get path to PFASS root directory
#PFASSTDIR="$(dirname $0)/.."

# Convert to absolute path
#PFASSTDIR="$(realpath $PFASSTDIR)"

#TESTDIR="$PFASSTDIR/test"

#echo "Using PFASSTDIR=$PFASSTDIR"
#echo "Using TESTDIR=$TESTDIR"


# Change to test directory
echo "Building libpfasst"
make clean
make -j || exit 1

echo "Building magpicard"
cd "test/magpicard"
make clean
make -j || exit 1

