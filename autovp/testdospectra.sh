#!/bin/bash

#gcc -o mkspecRAD -DDO9IONS mkspecRAD.c -lm

#gcc -o mkspecRAD -DDO9IONS -NONOISE mkspecRAD.c -lm

#autovpdir=/lustre/shuiyao/sci/autovp
sn=30
dvpix=6
dvres=0     
#ions='HI HeII CIII CIV OIV OVI NeVIII MgII SiIV'
ions='MgII'
taufact=0.62
#galmid=`awk '{print($5)}' <  ezw.MH$mass.zinbox.feb28b.all`
#taufile="specztau.test"
subfolder="test"
taudir=../$subfolder/
zquasar=0.26
ion=$1

#mkspecRAD $taufile $sn $zquasar $taufact $dvpix $dvres $ion >> auto.test.out
#mv $ion.raw $ion.$taufile

for ion in `echo $ions`
do
    for f in `ls $taudir | grep spec` 
    do
	echo mkspecRAD $taudir$f $sn $zquasar $taufact $dvpix $dvres $ion
	mkspecRAD $taudir$f $sn $zquasar $taufact $dvpix $dvres $ion >> out.mkspec.$subfolder

	echo mv $ion.raw $ion.cln
	mv $ion.raw $ion.cln

	echo autofit $ion
	autofit $ion >> out.autofit.$subfolder

	echo minfit $ion.pro
	minfit $ion.pro >> out.minfit.$subfolder

	echo mv $ion.vpm ../vpm/$subfolder/$ion/$f.vpm
	mv $ion.vpm ../vpm/$subfolder/$ion/$f.vpm
    done
done

