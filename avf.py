#!/usr/bin/env python

# **************************************************
# ..................................................
# Assemble variant files
# ..................................................
# **************************************************

# DESCRIPTION:
# The variant structure files are being assembled
# from sequence files. Each sequence file contains the
# data of one side chain and  every mutant has a
# fragment file.
# This script writes the BASH script which assembles
# the final mutant structure file(s).
# The variants are entered at the end of the script.

# CALLING SEQUENCE:
# $ python asv.py {1,3} -> seq.sh
# $ bash seq.sh

# PARAMETERS (search for by <<PARAM>>):
# - Reactant [= '1'] or product [= '3'] state files:

import sys
import fnmatch
import os


state=sys.argv[1] 
vars=sys.argv[2:]
if type(vars) == type(""):
    vars = [vars]
print sys.argv
print vars

# Available backbone sequence fragements.
chain = []
for i in range(1,501):
    for file in os.listdir("."):
        if fnmatch.fnmatch(file, 'seq-?%i-%s.pdb'%(i,state)):
            # Remove the 6 last chars (eg "-3.pdb")
            # And append to sequence of amino acids
            chain.append(file[0:-6])

def writeCatSeq(variant):
    # Convenience list of variants, replacing the '+'.
    varlist = variant.split('+')

    # Initialization of the bash script
    # Escaping BASH syntax.  
    # Resetting content of mutant structure file.
    numOfMutations = len(variant.split('+'))

    # <<PARAM>>:
    # - Reactant or product state file:        '?-res'
    # - Initial file or optimization job file: '-opt'.
    writeSeq  =\
            '# ' + '*'*50 + '\n' +\
            '# ' + '-'.join(variant.split('+')) + '\n' +\
            '# ' + str(numOfMutations) + '-fold mutant.\n' +\
            '# Resetting mutant file content.\n' +\
            'echo ' + '%'*50 + '\n'\
            'echo "Generating variant structure file of mutant:"\n' +\
            'echo ' + variant + '\n' +\
            'echo ' + '-'*50 + '\n'\
            'echo ' + '\n' +\
            'echo ' + '\n' +\
            'cat /dev/null > %s-res-' % state + '-'.join(varlist) + '-ste-ini.pdb\n' +\
            'for i in\\\n' +\
            '    {' 

    # Defining a list which contains only the numbers of the
    # residues that will be mutated.
    # First recast from 'G39A+L278A' form to ['G39A', 'L278A'],
    # then generate a list with only numbers ['39', '278'].
    varTmp = variant.split('+')
    varNum = [i[1:-1] for i in varTmp] 

    for s in chain:
        # Checking if the number of the residue is in the ones
        # to be mutated, if so, then the write sequence is
        # adjusted.
        if s[5:] in varNum:
            # Locate the index of the side chain which needs
            # to be mutated, then choose the corresponding index
            # in the varTmp list.
            ind=varNum.index(s[5:])
            writeSeq += 'frag-' + variant.split('+')[ind].lower() + ',\\\n'
        else:
            # Appending of the WT backbone.
            writeSeq += s + ',\\\n'

    # <<PARAM>>:
    # - Reactant or product state file:        '?-res'
    # - Initial file or optimization job file: '-opt'.
    writeSeq += '}\n' +\
            'do\n' +\
            '    cat $i-%s.pdb >> %s-res-' % (state, state) + '-'.join(varlist) + '-ste-ini.pdb\n' +\
            'done\n' +\
            'grep -v \'END\' ' + '%s-res-' % state + '-'.join(varlist) + '-ste-ini.pdb > tmp.pdb\n' +\
            'mv tmp.pdb ' + '%s-res-' % state + '-'.join(varlist) + '-ste-ini.pdb\n\n\n'
    return writeSeq

if __name__ == '__main__':
    #***********************
    # Reset the sequence script.
    seqFile = open('seq.sh', 'w')
    seqFile.close()
    seqFile = open('seq.sh', 'a') 
    #-----------------------


    #***********************
    # Generating the 'seq.sh' script
    for variant in vars:
        seqFile.write(writeCatSeq(variant))
    seqFile.close()
    #-----------------------

# EOF
#------------------------------------------------------
