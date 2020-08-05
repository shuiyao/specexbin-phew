import sys

errormsg = "Usage: "+sys.argv[0]+" modelname"
if(len(sys.argv) != 2):
    print errormsg
    sys.exit(1)
else:
    model = sys.argv[1]

# model = "l25n288-phew-m5"

basedir = "/proj/shuiyao/los/shortlos/"+model+"/"

mstrs = ["MH12", "MH11"]

#bpars = [10, 20, 30, 40, 50, 60, 80, 100, 150, 200, 250, 300]
bpars = [10, 50, 90, 130, 170, 210, 250, 290]
    
for mstr in mstrs:
    if(mstr == "MH11"): continue
    for bpar in bpars:
        bstr = str(bpar)+"kpc"
        datadir = basedir + mstr + "/" + bstr + "/"
        fnames = os.listdir(datadir)
        print datadir, ": ", len(fnames)

        for fname in fnames:
            fw = datadir + "specaimw." + fname[9:]
            fc = datadir + "specaimc." + fname[9:]
            fo = datadir + "specaimo." + fname[9:]

            fin1 = open(fw, "r")
            fin2 = open(fc, "r")
            fout = open(fo, "w")

            print "Writing:", fo
            for line in fin1:
                lineo = ""
                spt1 = line.split()
                spt2 = fin2.readline().split()
                for i in range(len(spt1)):
                    if(i in [7, 11, 15, 19, 23, 27, 31, 35, 39]):
                        lineo += "%5.3e " % ((float)(spt1[i]) + (float)(spt2[i]))
                        # lineo += "%5.3e " % ((float)(spt1[i]))
                    else:
                        lineo += spt1[i] + " "
                lineo = lineo + "\n"    # BE WARNED! THERE MUST BE A SPACE BEFORE THE \n. OTHERWISE SiIV will be 0
                fout.write(lineo)
            fin1.close()
            fin2.close()
            fout.close()
