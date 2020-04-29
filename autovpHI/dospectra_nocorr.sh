#!/bin/bash

gcc -o mkspecRAD -DDO9IONS mkspecRAD.c -lm

#gcc -o mkspecRAD -DDO9IONS -NONOISE mkspecRAD.c -lm

#autovpdir=/lustre/shuiyao/sci/autovp
sn=30
dvpix=6
dvres=0     
#ions='HI HeII CIII CIV OIV OVI NeVIII MgII SiIV'
ions='HI'
taufact=1.00
subfolder="z2tau1"
taudir=../z2/
zquasar=2.01

#mkspecRAD $taufile $sn $zquasar $taufact $dvpix $dvres $ion >> auto.test.out
#mv $ion.raw $ion.$taufile

mkdir ../vpm
mkdir ../vpm/$subfolder
mkdir ../raw
mkdir ../raw/$subfolder
for ion in `echo $ions`
do
    echo ================ $ion ================
    mkdir ../vpm/$subfolder/$ion
    mkdir ../raw/$subfolder/$ion
    for f in `ls $taudir | grep spec`
    do
	#mkspecRAD
	echo mkspecRAD $taudir$f $sn $zquasar $taufact $dvpix $dvres $ion
	./mkspecRAD $taudir$f $sn $zquasar $taufact $dvpix $dvres $ion >> out.mkspecRAD.$subfolder
	echo cp $ion.raw $ion.cln
	cp $ion.raw $ion.cln
	mv $ion.raw ../raw/$subfolder/$ion/$f.raw
	#autofit
	echo autofit $ion
	./autofit $ion >> out.autofit.$subfolder
	#minfit
	echo minfit $ion.pro
	./minfit $ion.pro >> out.minfit.$subfolder
	echo mv $ion.vpm ../vpm/$subfolder/$ion/$f.vpm
	mv $ion.vpm ../vpm/$subfolder/$ion/$f.vpm
    done
done

