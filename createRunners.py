import glob
import sys

year = sys.argv[1]

fileslists = glob.glob("filelists/{}/*.files".format(year))

filelists = [x.split('/')[2].split('.')[0] for x in fileslists]

dic = {}
for f in filelists:
    dic[f] = f+'.root'

f = open("runCount{}.sh".format(year),"w")
f.write("#!/bin/bash\n")
for k,v in dic.items():
    line = "root -l -b -q 'quickCount.C(\"{}\",\"{}\",{})'\n".format(k,v,year)
    f.write(line)
f.close()

f = open("runPUMC{}.sh".format(year),"w")
f.write("#!/bin/bash\n")
for k,v in dic.items():
    if 'SingleMuon' not in k:
        line = "root -l -b -q 'quickPUMC.C(\"{}\",\"{}\",{})'\n".format(k,v,year)
        f.write(line)
f.close()

f = open("runSelector{}.sh".format(year),"w")
f.write("#!/bin/bash\n")
for k,v in dic.items():
#     if year == "2016" and "M-105To160" in k:
#         continue
    line = "root -l -b -q 'quickScript.C(\"{}\",\"{}\",{})'\n".format(k,v,year)
    line2 = "mv {} ntupleFiles/{}/\n".format(v,year)
    f.write(line)
    f.write(line2)
f.close()