#!/bin/bash

model=$1

if [ -z $1 ]; then
    echo "Please specify modelname."
    echo "For example: bash dospectra.sh l25n288"
    exit -1
else
    fbase="/proj/shuiyao/los/"$model
    if [ -d $fbase ]; then
	echo "Working under: "$fbase
    else
	echo "Create workplace: "$fbase
	mkdir $fbase
    fi
    echo 
fi

# python add_clouds_to_spec_shortlos.py $model

# cp autovp.par.OVI ./src/autovp.par
cp autovp.par.HI ./src/autovp.par

sn=30
dvpix=6
dvres=0     
#ions='HI HeII CIII CIV OIV OVI NeVIII MgII SiIV'
# ions='CIV OVI NeVIII SiIV'
ions='HI'
taufact=0.31
subfolder="z0.5"
taudir=$fbase/$subfolder/
zquasar=0.51

# Re-compile mkspecRAD
pushd src/
gcc -o mkspecRAD -DDO9IONS mkspecRAD.c -lm
#gcc -o mkspecRAD -DDO9IONS -NONOISE mkspecRAD.c -lm

mkdir $fbase/vpm
mkdir $fbase/vpm/$subfolder
mkdir $fbase/raw
mkdir $fbase/raw/$subfolder

for ion in `echo $ions`
do
    echo ================ $ion ================
    mkdir $fbase/vpm/$subfolder/$ion
    mkdir $fbase/raw/$subfolder/$ion
    for f in `ls $taudir | grep specztauo`
    do
	#mkspecRAD
	echo mkspecRAD $taudir$f $sn $zquasar $taufact $dvpix $dvres $ion
	./mkspecRAD $taudir$f $sn $zquasar $taufact $dvpix $dvres $ion > out.mkspecRAD.$subfolder
	echo cp $ion.raw $ion.cln
	cp $ion.raw $ion.cln
	mv $ion.raw $fbase/raw/$subfolder/$ion/$f.raw
	#autofit
	echo autofit $ion
	./autofit $ion > out.autofit.$subfolder
	#minfit
	echo minfit $ion.pro
	./minfit $ion.pro > out.minfit.$subfolder
	echo mv $ion.vpm $fbase/vpm/$subfolder/$ion/$f.vpm
	mv $ion.vpm $fbase/vpm/$subfolder/$ion/$f.vpm
    done
done

popd
