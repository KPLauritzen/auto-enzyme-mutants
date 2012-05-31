for mutant in $(cat mutantList.dat)
do
    ROT=1
    while [ -e 3-$mutant-rot$ROT.pdb.out ]
    do
        grep -i "final heat" 3-$mutant-rot$ROT.pdb.out |
          sed 's:.* = ::' | sed 's:KJ/MOL::' >> $mutant.dat
        let ROT=ROT+1
    done
    lowest=$(cat $mutant.dat | sort | tail -1)
    bestRot=`grep -n -e $lowest $mutant.dat | sed 's/\(.*\):.*/\1/'`
    mutLower=`echo $mutant | awk '{print tolower($0)}'`
    cp frag-$mutLower-rot$bestRot-3.pdb frag-$mutLower-3.pdb
    cp frag-$mutLower-rot$bestRot-3.pdb frag-$mutLower-1.pdb
    rm $mutant.dat
    python avf.py 1 $mutant
    bash seq.sh
    python avf.py 3 $mutant
    bash seq.sh
done
