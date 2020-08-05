import sys

errormsg = "Usage: "+sys.argv[0]+" modelname"
if(len(sys.argv) != 2):
    print errormsg
    sys.exit(1)
else:
    model = sys.argv[1]

angles = range(10, 81, 1)

for ai in range(10, 81, 1):
# for ai in [13]:    
    fw = "/proj/shuiyao/los/"+model+"/z0.5/specztauw."+model+"."+str(ai)+".50_0"
    fc = "/proj/shuiyao/los/"+model+"/z0.5/specztauc."+model+"."+str(ai)+".50_0"
    fo = "/proj/shuiyao/los/"+model+"/z0.5/specztauo."+model+"."+str(ai)+".50_0"

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
