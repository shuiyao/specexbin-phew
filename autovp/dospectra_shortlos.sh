#!/bin/bash

model=$1

if [ -z $1 ]; then
    echo "Please specify modelname."
    echo "For example: bash dospectra.sh l25n288"
    exit -1
else
    fbase="/proj/shuiyao/los/shortlos/"$model
    if [ -d $fbase ]; then
	echo "Working under: "$fbase
    else
	echo "Create workplace: "$fbase
	mkdir $fbase
    fi
    echo 
fi

python add_clouds_to_spec.py $model

cp autovp.par.OVI ./src/autovp.par
# cp autovp.par.HI ./src/autovp.par

sn=30
dvpix=6
dvres=0     
#ions='HI HeII CIII CIV OIV OVI NeVIII MgII SiIV'
#blist='10 20 30 40 50 60 80 100 150 200 250 300'
blist='10 50 90 130 170 210 250 290'
ions='OVI MgII CIV'
taufact=1.0
suffix="z0.25"
zquasar=0.26

# Re-compile mkspecRAD
pushd src/
gcc -o mkspecRAD -DDO9IONS mkspecRAD.c -lm
#gcc -o mkspecRAD -DDO9IONS -NONOISE mkspecRAD.c -lm

mkdir $fbase/raw
mkdir $fbase/vpm

for ion in `echo $ions`
do
    echo ================ $ion ================
    mkdir $fbase/raw/$ion
    mkdir $fbase/vpm/$ion    
#     for mh in MH11 MH12
    for mh in MH12
    do
	mkdir $fbase/raw/$ion/$mh
	mkdir $fbase/vpm/$ion/$mh	
	for b in `echo $blist`
	do
    	    taudir=$mh/$b"kpc"
	    mkdir $fbase/raw/$ion/$taudir
	    mkdir $fbase/vpm/$ion/$taudir	    
	    for f in `ls $fbase/$taudir | grep specaimo`
	    do
		#mkspecRAD
 		echo mkspecRAD $fbase/$taudir/$f $sn $zquasar $taufact $dvpix $dvres $ion
 		./mkspecRAD $fbase/$taudir/$f $sn $zquasar $taufact $dvpix $dvres $ion > out.mkspecRAD.$suffix
 		cp $ion.raw $fbase/raw/$ion/$taudir/$f.raw
 		echo mv $ion.raw $ion.cln
 		mv $ion.raw $ion.cln
		#autofit
		echo autofit $ion
		./autofit $ion > out.autofit.$suffix
		#minfit
		echo minfit $ion.pro
		./minfit $ion.pro > out.minfit.$suffix
		echo mv $ion.vpm $fbase/vpm/$ion/$taudir/$f.vpm
		mv $ion.vpm $fbase/vpm/$ion/$taudir/$f.vpm
	    done
    	done
    done
done

popd
