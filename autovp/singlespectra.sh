#!/bin/bash

rm mkspecRAD
#gcc -o mkspecRAD -DDO9IONS mkspecRAD.c -lm

gcc -o mkspecRAD -DDO9IONS -DNONOISE mkspecRAD.c -lm

#autovpdir=/lustre/shuiyao/sci/autovp
sn=30
dvpix=6
dvres=0     
#ions='HI HeII CIII CIV OIV OVI NeVIII MgII SiIV'
ions='HI OVI MgII CIV'
#ions='HI'
taufact=1.0
subfolder="z0.25"
zquasar=0.26

# f="../ztauw.37.5_0"
# prefix="pp"
# f="../ztau.37.5_0"
# prefix="sm"
#f="../shortlos.p50n288o5pp.10kpc"
#prefix="pp"
f="../specaimw.shortlos.h727.pp"
prefix="pp"
# f="../specaimw.shortlos.h727.ppeasy"
# prefix="ppe"



for ion in `echo $ions`
do
    echo ================ $ion ================
    #mkspecRad
    ./mkspecRAD ./$f $sn $zquasar $taufact $dvpix $dvres $ion >> out.mkspecRAD.$subfolder
    echo cp $ion.raw $ion.cln
    cp $ion.raw $ion.cln
    cp $ion.cln $prefix.$ion.cln
    #autofit
    echo autofit $ion
    ./autofit $ion >> out.autofit.$subfolder
    #minfit
    echo minfit $ion.pro
    ./minfit $ion.pro >> out.minfit.$subfolder
done    

