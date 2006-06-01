#!/usr/bin/env python

import synthlisa
import Numeric                                
import pyx

import sys
import os
import re

def multiplot(spectra,titles,filename,number,loglog=False):
    # will need more styles!
    style = [pyx.graph.style.line([pyx.color.rgb.red]),
             pyx.graph.style.line([pyx.color.rgb.green]),
             pyx.graph.style.line([pyx.color.rgb.blue]),
             pyx.graph.style.line([pyx.color.rgb(1.0,1.0,0.0)]),
             pyx.graph.style.line([pyx.color.rgb(1.0,0.0,1.0)]),
             pyx.graph.style.line([pyx.color.rgb(0.0,1.0,1.0)])]

    if loglog:
        g = pyx.graph.graphxy(width=16,
                              key=pyx.graph.key.key(pos="tl"),
                              x=pyx.graph.axis.log(),
                              y=pyx.graph.axis.log())
    else:
        g = pyx.graph.graphxy(width=16,
                              key=pyx.graph.key.key(pos="tl"))

    # if type(spectra) == list:
    #     for i in range(0,len(spectra)):
    #         g.plot(pyx.graph.data.list(map(tuple,spectra[i]),
    #                                    x=1,y=2,title=titles[i]),
    #                [style[i]])
    # elif type(spectra) == Numeric.ArrayType:
    
    for i in range(1,Numeric.shape(spectra)[1]):
        g.plot(pyx.graph.data.list(map(tuple,spectra),
                                   x=1,y=i+1,title=titles[i-1]),
               [style[i-1]])

    if type(filename) == list:
        try:
            thefilename = filename[number]
        except:
            print "multiplot: need more output files!"
    else:
        if '.eps' in filename:
            thefilename = re.sub('\.eps$','-' + str(number) + '.eps',filename)
        elif '.pdf' in filename:
            thefilename = re.sub('\.pdf$','-' + str(number) + '.pdf',filename)
        else:
            thefilename = filename + '-' + str(number) + '.eps'
            
    if '.pdf' in thefilename:
        g.writePDFfile(thefilename)
    else:
        g.writeEPSfile(thefilename)

    return thefilename


r = synthlisa.readXML(sys.argv[1])

spectra = r.getTDIFrequencySeries()
timeseries = r.getTDITimeSeries()

# number 

if len(sys.argv) == 2:
    basename = re.sub('\.xml$','',sys.argv[1])
else:
    basename = sys.argv[2:]

counter = 0
files = []

for spectrum in spectra:
    if counter < len(basename):
        files.append(multiplot(spectrum['Data'],spectrum['Vars'][1:],basename,counter,loglog=True))
        counter = counter + 1

for series in timeseries:
    if counter < len(basename):
        files.append(multiplot(series['Data'],series['Vars'][1:],basename,counter,loglog=False))
        counter = counter + 1

if sys.platform == 'darwin':
    for file in files:
        os.system('open ' + file)
elif 'linux' in sys.platform:
    for file in files:
        os.system('gv ' + file + ' &')

# do some magic to do multiple plots by collating single spectra,
# or by using one multiple spectrum

# collect = []
# titles = []

# for spectrum in spectra:
#     if Numeric.shape(spectrum['Data'])[1] > 2:
#         if collect:
#             multiplot(collect,titles,basename,counter)
#             counter = counter + 1
#             collect = []
#             titles = []
#         multiplot(spectrum['Data'],spectrum['Vars'][1:],basename,counter)
#         counter = counter + 1
#     else:
#         collect.append(spectrum['Data'])
#         titles.append(spectrum['Vars'][1])
#
# if collect:
#     multiplot(collect,titles,basename,counter)
    
# ??? encapsulate in lisautils.py?
# ??? legend?
