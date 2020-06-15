#!/bin/bash
#parallel contspecexbin_v8 {1} {2} {3} {4} {5} {6} ::: snap_p50n288ezw15 ::: 0.0 ::: 6.0 ::: 50.0 ::: 1.0 :::: angles_10_24

modelname=l50n288-phewoff
zbeg=0.0
zend=0.5
lbox=50.0
ftau=1.0
mc=2.e37 # PhEW Cloud Mass

cd .. 
while read line 
do
    echo "./contspecexbin "$modelname $zbeg $zend $lbox $ftau $line
    # gdb --args ./contspecexbin $modelname $zbeg $zend $lbox $ftau $line
    ./contspecexbin $modelname $zbeg $zend $lbox $ftau $line $mc
done < "angles.dat"
