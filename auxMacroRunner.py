#!/usr/bin/env/python

from subprocess import Popen

#see inside outputFiles for all the possible configs
histos = ['p','pt','eta','phi']
scenarios = ['0010','0100','1000']
regions = ['barrel24d','endcap','overlap']
dirs = ['GLB','FMS']

for scenario in scenarios:
    for dir in dirs:
        for region in regions:
            for hist in histos:
                Popen('root -b -q -l \'superimpose.C(\"'+hist+'\",\"'+scenario+'\",\"'+region+'\",\"'+dir+'\")\'',shell=True).wait()
