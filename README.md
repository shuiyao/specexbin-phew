# SPECEXBIN:

Generate mock line of sight (LoS) quasar absorption spectra from cosmological simulations. It is made compatible with the GIZMO code + PhEW model. It contains AUTOVP, a software that analyzes the spectra data and perform profile fitting.

## Files and Parameters
================================
Here is a list of user defined files and parameters that you might need to reconfigure on a different machine.

_./specions_i9_gizmo.dat_: The ions data, including atomic weight, solar mass fraction, corresponding column in the P.Metallicity field in the simulation.
_./Makefile_:
  - SHORTSPEC: Turn on to enable shortLoS, otherwise generate random LoS
_./specexbindefs.h_:
  - FOLDER_IONS: The folder that contains data tables that are necessary
  - FOLDER_OUTPUT: The default output folder. Since thousands of files will be dumped into this folder. We better choose a "dump site" that is exclusively used for temporarily storing the data.
_./SCRIPT/angles_: The angle files, used for generating random LoS
_./SCRIPT/tabs_: The files that point to the locations of the snapshots from the simulation.
_./SCRIPT/LoS_: The targeted LoS (shortLoS) that was generated.

## Generate Random LoS with SPECEXBIN
================================
Suppose now we want to generate random LoS for a simulation named $modelname.

1. Preparation
  - In SCRIPT/tabs/, create the shortlos_$modelname.tab file. You can edit from the sample.tab file.
  - Check SHORTSPEC, OUTPUT_LOCAL_FOLDER and specexbindefs.h
  - Generate the .slm files and .sh files for a given model. For example:
  ```
  python create_batch_files.py l25n288 2.0e38
  ```
  - Check the .slm and .sh file, then
2. Run SPECEXBIN with
```
bash batch.sh
```
3. Distribute Output Files
  - The output files are temporarily stored at FOLDER_OUTPUT
  - 


4. (PhEW only) Add Clouds to the Spectra

In the autovp folder, use
```
python add_clouds_to_spec.py $modelname
```

6. Do AutoVP

```
bash dospectra.sh
```

7. Compile VPM files
  - spec_compile.py and spec_compile_HI.py

## Generate ShortLOS with SPECEXBIN.
================================
1. Generate LoS info
  - mklosfile.py
    - Purpose: Create LoS data in LoS/ [A]
    - Parameters: modelname, MASSBIN=12, NHALO=250
  - In tabs/, create the shortlos_$modelname.tab file [B]
2. Specexbin
  - Compile with SHORTSPEC, check OUTPUT_LOCAL_FOLDER and specexbindefs.h
  - shortlos.slm -> batch_shortlos.sh
    - Specify modelname, mcinit, lbox. Takes [A], [B] as input
    - By default, output in /proj/shuiyao/los/shortlos/temp
3. Distribute specaim files
  - distributeLOS.sh
    - Calls createlosfolder.sh
    - Calls distributeLOS.py
      - Create a $modelname.page file [C]
    - move files from /proj/shuiyao/los/shortlos/temp to distinations according to [C]
4. Copy autovp/ and autovpHI/ to $modelname folder
5. (PhEW only) Add clouds to the spectra
  - autovp/add_clouds_to_spec.py
6. AutoVP
  - Use dospectra.sh
7. Compile VPM files
  - spec_compile.py and spec_compile_HI.py
