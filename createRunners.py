import glob
import sys
import os

year = sys.argv[1]

fileslists = glob.glob("filelists/{}/*.files".format(year))

filelists = [x.split('/')[2].split('.')[0] for x in fileslists]

dic = {}
for f in filelists:
    dic[f] = f+'.root'

f = open("runCount{}.sh".format(year),"w")
f.write("#!/bin/bash\n")
for k,v in dic.items():
    if "SingleMuon" not in k:
        line = "root -l -b -q 'quickCount.C(\"{}\",\"{}\",{})'\n".format(k,v,year)
        f.write(line)
f.close()
os.system("chmod 755 {}".format("runCount{}.sh".format(year)))

f = open("runCorr{}.sh".format(year),"w")
f.write("#!/bin/bash\n")
for k,v in dic.items():
    if "SingleMuon" not in k:
        line = "root -l -b -q 'quickCorr.C(\"{}\",\"{}\",{})'\n".format(k,v,year)
        f.write(line)
#f.write("hadd resources/data/allData{}.root resources/corrections/{}/SingleMuon*.root\n".format(year,year))
f.close()
os.system("chmod 755 {}".format("runCorr{}.sh".format(year)))

f = open("runSelector{}.sh".format(year),"w")
f.write("#!/bin/bash\n")
for k,v in dic.items():
    line = "root -l -b -q 'quickScript.C(\"{}\",\"{}\",{})'\n".format(k,v,year)
    line2 = "mv {} ntupleFiles/{}/\n".format(v,year)
    f.write(line)
    f.write(line2)
f.close()
os.system("chmod 755 {}".format("runSelector{}.sh".format(year)))

f = open("runHistos{}.sh".format(year),"w")
f.write("#!/bin/bash\n")
for k,v in dic.items():
    line = "root -l -b -q 'histoScript.C(\"{}\",{})'\n".format(v,year)
    f.write(line)
f.write("mkdir histoFiles/{}/data\n".format(year))
f.write("mv histoFiles/{}/SingleMuon*.root histoFiles/{}/data/\n".format(year,year))
f.write("hadd histoFiles/{}/allData{}.root histoFiles/{}/data/*.root\n".format(year,year,year))
f.close()
os.system("chmod 755 {}".format("runHistos{}.sh".format(year)))
