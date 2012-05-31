import __main__
__main__.pymol_argv = ['pymol','-qc'] # Pymol: quiet and no GUI
from time import sleep
import pymol
pymol.finish_launching()
from pymol import cmd
from pymol import stored
from pymol.exporting import _resn_to_aa as one_letter
import os
from os.path import splitext
import sys
import fnmatch

#PARAMETERS
state = "3"

# Filenames for initial and end-state of the reaction
obj3 = '3-wt-opt.pdb'
obj1= '1-wt-opt.pdb'

# For convience. Allows for defining mutations like: '140' : allAminoAcids
allAminoAcids = [ 'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN',
                  'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE',
                  'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL' ]
# Define which mutation should be done. e.g. '140': ['ARG', 'LYS']
mutations = {
'140': ['ARG', 'LYS']
}

# commandname for MOPAC submit-script
mopacCommand = "echo"

# *****************************************************************
def seq(state, selection="name ca or resn hoh or resn lig"):
    print "Generating seqs."
    cmd.select("prot", selection)
    while cmd.pop("_tmp", "prot"):
        cmd.iterate("_tmp", "stored.x=(resn,resv)")
        #print stored.x[0], stored.x[1]
        # Special case 1: Waters.
        if stored.x[0] == 'HOH':
            filename = 'seq-x%s-%s.pdb' % (stored.x[1], state)
        # Special case 2: Substrate.
        elif stored.x[0] == 'LIG':
            filename = 'seq-x%s-%s.pdb' % (stored.x[1], state)
        # Other: protein back-bone.
        else:
            filename = 'seq-%s%d-%s.pdb' \
            % (one_letter[stored.x[0]].lower(), stored.x[1], state)
        cmd.save(filename, "byres _tmp")
    cmd.delete('_tmp prot')



# ************************************************************
def setup(obj):

    # Various set up
    pwd = os.getcwd() 
    #print "os.getcwd()", os.getcwd()
    cmd.do('wizard mutagenesis')
    cmd.do('refresh_wizard')

    # Save residue names and numbers.
    orig_sequence = setNames(obj)
    return pwd, orig_sequence
# ------------------------------------------------------------


# ************************************************************
# 'state=state': The first variable is the variable used within
# the scope of this function. The second variable is the one
# in the global scoped and defined at the top of the module.
def frag(state=state, obj=obj3): 
    pwd, orig_sequence = setup(obj)
    stored.rotamerDict = {}
    # Add and retain hydrogens
    cmd.get_wizard().set_hyd("keep") 

    # Run over all sites where to mutate
    for site in mutations.keys():

        variants = mutations[site]

        # Run over all variants.
        for variant in variants:
            cmd.load(obj) 

            cmd.do('wizard mutagenesis')
            cmd.do('refresh_wizard')
            cmd.get_wizard().do_select("(%s/)" % site)
            cmd.do("cmd.get_wizard().set_mode('%s')"%variant)

            # Get the number of available rotamers at that site
            # Introduce a condition here to check if 
            # rotamers are requested. 
            # <<OPTION>>
            print variant, "variant"
            nRots = getRots(site, variant)
            nRots = 2
            stored.rotamerDict[str(site)+getOne(variant)] = nRots

            cmd.rewind()
            for i in range(1, nRots + 1): 
                
                cmd.get_wizard().do_select("(" + site + "/)")
                cmd.frame(i)
                cmd.get_wizard().apply()

                # Optimize the mutated sidechain
                #<<OPTION>>
                print "Sculpting."
                localSculpt(obj, site)

                # Protonation of the N.
                cmd.do("select n%d, name n and %d/" % (int(site), int(site)))
                cmd.edit("n%d" % int(site), None, None, None, pkresi=0, pkbond=0)
                cmd.do("h_fill")

                # Protonation of the C.
                cmd.do("select c%d, name c and %d/" % (int(site), int(site)))
                cmd.edit("c%d" % int(site), None, None, None, pkresi=0, pkbond=0)
                cmd.do("h_fill") 

                # Definition of saveString
                saveString  = '%s/' % pwd
                saveString += 'frag-' + getOne(orig_sequence[site]).lower() +\
                               site + getOne(variant).lower() + '-rot%i-%s.pdb, ' \
                               % (i,state) +'((%s/))' % site
                #print saveString 
                cmd.do('save %s' % saveString)
            cmd.do('delete all') 
            cmd.set_wizard('done')
    print "Frag is all done"

# ------------------------------------------------------------


# ************************************************************
# Convenience Functions
def getRots(site, variant): 
    cmd.get_wizard().set_mode("\""+variant+"\"")
    
    # Key lines 
    # I dont know how they work, but they make it possible.
    # Jason wrote this: If you just write "site" instead of
    #                   "(site)", PyMOL will delete your
    #                   residue. "(site)" makes it an
    #                   anonymous selection.
    #print 'getRots'
    cmd.get_wizard().do_select("(" + str(site) + "/)")
    nRot = cmd.count_states("mutation") 
    return nRot 

def setNames(obj):
    orig_sequence = {}
    cmd.load(obj) 
    cmd.select("prot", "name ca")
    cmd.do("stored.names = []")
    cmd.do("iterate (prot), stored.names.append((resi, resn))")
    for i in stored.names:
        orig_sequence[i[0]] = i[1] 
    cmd.do('delete all') 
    stored.orig_sequence = orig_sequence
    return orig_sequence



# Credit: Thomas Holder, MPI
# CONSTRUCT: - 'res'
#            - 'cpy'
#            -
def localSculpt(obj, site):
    res = str(site)
    cmd.protect('(not %s/) or name CA+C+N+O+OXT' % (res))
    print "Activating Sculpting."
    cmd.sculpt_activate(obj[:-4]) 
    cmd.sculpt_iterate(obj[:-4], cycles=5000) 
    cmd.sculpt_deactivate(obj[:-4])
    cmd.deprotect() 


def getOne(three):
    trans = { 
       'ALA':'A',
       'ARG':'R',
       'ASN':'N',
       'ASP':'D',
       'CYS':'C',
       'GLU':'E',
       'GLN':'Q',
       'GLY':'G',
       'HIS':'H',
       'ILE':'I',
       'LEU':'L',
       'LYS':'K',
       'MET':'M',
       'PHE':'F',
       'PRO':'P',
       'SER':'S',
       'THR':'T',
       'TRP':'W',
       'TYR':'Y',
       'VAL':'V'
       } 
    return trans[three]
# ------------------------------------------------------------


# ************************************************************
# Expose to the PyMOL shell
cmd.extend('setup', setup)
cmd.extend('frag', frag)
cmd.extend('getRots', getRots)
cmd.extend('localSculpt', localSculpt)
cmd.extend('seq', seq)
# ------------------------------------------------------------

######################
# Modified avf.py
#####################
            
def writeCatSeq(variant):
    # Convenience list of variants, replacing the '+'.
    varlist = variant.split('+')

    # Initialization of the bash script
    # Escaping BASH syntax.  
    # Resetting content of mutant structure file.
    numOfMutations = len(variant.split('+'))

    # Defining a list which contains only the numbers of the
    # residues that will be mutated.
    varTmp = variant.split('+')
    varNum = [i[1:-1] for i in varTmp]
    site = varNum[0]
    nRots = stored.rotamerDict[variant[1:]]
    
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
            'for ROT in $(seq %s)'%nRots +'\n' +\
            'do \n' +\
            'cat /dev/null > %s-' % state + '-'.join(varlist) + '-rot$ROT.pdb\n' +\
            'echo "1scf mozyme cutoff=15 \n\n" >' +\
            '%s-' % state + '-'.join(varlist) + '-rot$ROT.pdb \n' +\
            'for i in\\\n' +\
            '    {' 


    for s in chain:
        # Checking if the number of the residue is in the ones
        # to be mutated, if so, then the write sequence is
        # adjusted.
        if s[5:] in varNum:
            # Locate the index of the side chain which needs
            # to be mutated, then choose the corresponding index
            # in the varTmp list.
            ind=varNum.index(s[5:])
            writeSeq += 'frag-' + variant.split('+')[ind].lower() +"-rot$ROT"+ ',\\\n'
        else:
            # Appending of the WT backbone.
            writeSeq += s + ',\\\n'

    # <<PARAM>>:
    # - Reactant or product state file:        '?-res'
    # - Initial file or optimization job file: '-opt'.
    writeSeq += '}\n' +\
            'do\n' +\
            '    cat $i-%s.pdb >> %s-' % (state, state) + '-'.join(varlist) +\
                '-rot$ROT.pdb\n' +\
            'done\n' +\
            'grep -v \'END\' ' + '%s-' % state + '-'.join(varlist) +\
                '-rot$ROT.pdb > tmp.pdb\n' +\
            'mv tmp.pdb ' + '%s-' % state + '-'.join(varlist) + '-rot$ROT.pdb\n' +\
            '%s %s-' %(mopacCommand, state) + '-'.join(varlist) + '-rot$ROT.pdb\n' +\
            'done\n' +\
            'echo ' + '-'.join(varlist) + ' >> mutantList.dat\n'
    return writeSeq

if __name__ == '__main__':
    pymol.cmd.reinitialize()
    pymol.cmd.load(obj1)
    pymol.cmd.do("run create-mutant-fragments.py")
    # Sequence the 1-state
    pymol.cmd.do("seq 1")
    sleep(2)
    pymol.cmd.do("delete all")
    # Sequence the 3-state
    pymol.cmd.load(obj3)
    pymol.cmd.do("seq 3")
    sleep(2)
    pymol.cmd.sync()
    pymol.cmd.refresh()
    # Create the mutant fragments
    pymol.cmd.do("frag")
    pymol.cmd.refresh()
    pymol.cmd.sync(100000.0 * len(mutations))

    # This assumes that seq_files.py has already been run.
    # Define the backbone chain
    chain = []
    for i in range(1,551):
        for file in os.listdir("."):
            if fnmatch.fnmatch(file, 'seq-?%i-%s.pdb'%(i,state)):
                # Remove the 6 last chars (eg "-3.pdb")
                # And append to sequence of amino acids
                chain.append(file[0:-6])
    
    #***********************
    # Reset the sequence script.
    seqFile = open('seq.sh', 'w')
    seqFile.close()
    seqFile = open('seq.sh', 'a') 
    #-----------------------

    # Get a list of the mutations in one-letter, site, one-letter format. E.g. L140K
    vars = []
    for site,mutationList in mutations.items():
        for mutant in mutationList:
            vars.append(getOne(stored.orig_sequence[site]) + str(site) + getOne(mutant))
    print vars
    #***********************
    # Generating the 'seq.sh' script
    for variant in vars:
        seqFile.write(writeCatSeq(variant))
    seqFile.close()
    #-----------------------
    os.system("bash seq.sh")
    pymol.cmd.do("quit")
# EOF
#------------------------------------------------------
