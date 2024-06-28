# encoding=utf8
import sys

''' This script imports a list of molecules, sort them by alphabet order,
and generate the code section for GOBASIC catalog paths with the correct
Matlab format. It also generates the corresponding label of the fitting
components for GOBASIC output files. '''

listname = sys.argv[1]

with open(listname+'.list', 'r') as complist:
    compnames = complist.readlines()
compnames.sort()

with open(listname+'.m', 'w+') as output:
    output.write('molpath = {\n')
    for name in compnames:
        output.write('           \'../CatalogFile/{0}.txt\';\n'.format(name.strip('\n')))
    output.write('           };\n\n')
    output.write('parpath = {\n')
    for name in compnames:
        output.write('           \'../CatalogFile/{0}par.par\';\n'.format(name.strip('\n')))
    output.write('           };\n\n')
    output.write('parco = {\n')
    for name in compnames:
        output.write('           \'../CatalogFile/{0}par.co\';\n'.format(name.strip('\n')))
    output.write('           };\n\n')
    output.write('label = {\n')
    for name in compnames:
        output.write('         \'{}\';\n'.format(name.strip('\n')))
    output.write('         };')
