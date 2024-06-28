# encoding = utf8

''' This script is used for picking out spectra of molecular components in
a GOBASIC output file. The user inputs the names of molecules(tags) to be
picked. This script searches for these components in the header information,
and picks them out to generate a new output csv file. '''
 
import numpy as np
import sys
import argparse

# parse arguments
parser = argparse.ArgumentParser(description=__doc__,
                                 epilog='--- Luyao Zou, June 2015 ---')
parser.add_argument('file', nargs=1, help='input file name')
parser.add_argument('-o', '--out', nargs=1, help='output file name')
args = parser.parse_args()

in_name = args.file[0]
if args.out:
    out_name = args.out[0]
else:
    out_name = 'GOBASIC_Picked.csv'

mol_name = input('Picking molecule names (separate by space): ')

# read headers
with open(in_name, 'r') as fit_doc:
    version = fit_doc.readline()
    header = fit_doc.readline()

# split header and molecule names into lists
col_name_list = tuple(header.split(','))
mol_name_list = mol_name.split()
cols = []

# loop through typed-in molecule names and match column indices in the
# fit file
for mol in mol_name_list:
    # check duplication
    count = col_name_list.count(mol)
    # if name matches, record the column number
    if count:
        for n in range(count):
            cols.append(col_name_list.index(mol)+n)
cols.sort()

# print some information on screen
print('>>>>>>>>>>>>>>>>>>')
print('>> Find {:d} components.'.format(len(cols)))
if not cols:
    exit()

sys.stdout.write('>> Processing data >')
sys.stdout.flush()

# load data
fit = np.loadtxt(in_name, comments='#', delimiter=',', skiprows=2)
# prepare for output data table
out_table = np.empty([len(fit),len(cols)+4])
# copy obsfreq, obsint and previousfit into output data table
out_table[:,(0,1,2)] = np.compress([True, True, True], fit, axis=1)
out_header = (version.strip('\n') +
              ', and picked by python3.4 script PickGOBASIC.py\n'
              + 'obsfreq,obsint,previousfit,')
# copy the corresponding data columns in the fit file to output table
for n in cols:
    # pick out data column
    sys.stdout.write('>')
    sys.stdout.flush()
    out_table[:,cols.index(n)+3] = fit[:,n]
    out_header += col_name_list[n] + ','
# calculate sum_fit
out_table[:,-1] = np.sum(out_table[:,3:-1], axis=1)

# set output file headers
out_header += 'sum_fits'

# output
np.savetxt(out_name, out_table, fmt='%.10g', newline='\n', delimiter=',',
           header=out_header, comments='')
print('\n>> Output file saved as {:s}'.format(out_name))

