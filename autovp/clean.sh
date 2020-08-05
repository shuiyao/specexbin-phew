#!/bin/bash

# Clean out the temporary files.
pushd src
echo "Removing Files in: "
pwd
rm -f ./*.ewlim
rm -f ./*.cln
rm -f ./*.pro
rm -f ./*.min
rm -f ./*.vpm
rm -f ./*.raw
rm -f ./*.mod

rm -f ./out.*
popd
