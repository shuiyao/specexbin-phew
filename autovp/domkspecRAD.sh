#!/bin/bash

#gcc -o mkspecRAD -DDO9IONS mkspecRAD.c -lm

gcc -o mkspecRAD -DDO9IONS -DNONOISE mkspecRAD.c -lm

#autovpdir=/lustre/shuiyao/sci/autovp
sn=30
dvpix=6
dvres=0     
#ions='HI HeII CIII CIV OIV OVI NeVIII MgII SiIV'
#ions='HI CIV OVI MgII'
ions='HI'
taufact=0.62
#galmid=`awk '{print($5)}' <  ezw.MH$mass.zinbox.feb28b.all`
#taufile="specztau.test"
subfolder="test"
taudir=../$subfolder/
zquasar=2.01
ion=$1

#mkspecRAD $taufile $sn $zquasar $taufact $dvpix $dvres $ion >> auto.test.out
#mv $ion.raw $ion.$taufile

for ion in `echo $ions`
do
    echo ================ $ion ================
    for f in `ls $taudir | grep specztau`
    do
	#mkspecRAD
	echo mkspecRAD $taudir$f $sn $zquasar $taufact $dvpix $dvres $ion
	mkspecRAD $taudir$f $sn $zquasar $taufact $dvpix $dvres $ion >> out.mkspecRAD.$subfolder
# 	echo mv $ion.raw $ion.cln
# 	mv $ion.raw $ion.cln
# 	#autofit
# 	echo autofit $ion
# 	autofit $ion >> out.autofit.$subfolder
# 	#minfit
# 	echo minfit $ion.pro
# 	minfit $ion.pro >> out.minfit.$subfolder
# 	echo mv $ion.vpm ../vpm/$subfolder/$ion/$f.vpm
# 	mv $ion.vpm ../vpm/$subfolder/$ion/$f.vpm
    done
done

