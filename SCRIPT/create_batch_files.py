finname = "sample.slm"
fshname = "sample.sh"
foutname = "test.slm"

model = "l25n288-phew-m5"
mcinit = "2.0e38"

angles1 = [10., 20., 30., 40., 50., 60., 70.]
angles2 = [20., 30., 40., 50., 60., 70., 80.]
for fi in range(len(angles1)):
    a1, a2 = angles1[fi], angles2[fi]
    fin = open(finname, "r")
    foutname = model + "_" + str(a1)[:2] + "-" + str(a2)[:2] + ".slm"
    anglename = "angles_" + str(a1)[:2] + "_" + str(a2)[:2] + ".dat"
    fout = open(foutname, "w")
    print "Write: ", foutname
    for i in range(22):
        line = fin.readline()
        fout.write(line)
    line = "bash " + model + "_" + str(a1)[:2] + "-" + str(a2)[:2] + ".sh"
    fout.write(line)
    fin.close()
    fout.close()

    fin = open(fshname, "r")
    foutname = model + "_" + str(a1)[:2] + "-" + str(a2)[:2] + ".sh"
    fout = open(foutname, "w")
    print "Write: ", foutname    
    for i in range(3):
        line = fin.readline()
        fout.write(line)
    line = "modelname=" + model + "\n"
    fout.write(line)
    line = "mcinit=" + mcinit + "\n"
    fout.write(line)    
    fin.readline() # skip the modelname= line
    fin.readline() # skip the mcinit= line
    for i in range(11):
        line = fin.readline()
        fout.write(line)
    line = "done < " + anglename
    fout.write(line)
    fin.close()
    fout.close()
        
