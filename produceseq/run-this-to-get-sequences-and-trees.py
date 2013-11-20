import time

from generate_fasta_newick import *
#from runrealecotypes import *
from produce_sequences import *
#from runecosim import *
#from runrandom import *
#from rungmyc import * #just the little nelder-mead issue
#from runbaps import *
#from runaml import *

#from backtofasta import backtofasta

import sys, os
import os.path

runflag='A'

start= time.time()



likelihoods = [0]
path_to_fastTree = "./setup/FastTree"
#Test variables
drift,lengthseq = '1.0e25', '1500'
init_hab_num = 2
genes = 7
omega = .19
sigma = 1.1
loops = 1
#nus=[5-10-15-20]



# Aaron: CHANGE INPUT FROM HERE:


#YOU NEED TO CHANGE PATH IN 2 DIFFERENT SPOT
# 1: CHANGE THE PATH FOR QUICK_GENERATE BELOW
# 2: CHANGE THE PATH IN PRODUCESEQUENCES.PY


nus=['3'] # NUMBER OF SEQUENCES
npops = [1] # NUMBER OF ECOTYPES IN THERE


totalrun= len(nus)*len(npops)*loops

ct=1

print 'Total '+str(totalrun)+' runs'
print ''

for i in range(0,loops):
    print 'Loop '+str(ct)
    for nu in nus:
        print 'Test parameters:'
        print 'nu='+nu
        for npop in npops:
            print 'npop='+str(npop)
            print ''
            print 'Producing test sequence..'
            produce_sequences(omega, sigma, npop, drift, nu, lengthseq)
            print ''
            print 'Generating newick tree..'

            # CHANGE THE PATH HERE TO YOUR PRODUCESEQ FOLDER
            os.chdir('/Users/Wei/Desktop/produceseq/')
            quick_generate('./producesequences/produceseqs.dat', likelihoods[0], path_to_fastTree)
            print ''
            id='../output/'+runflag+':'+nu+':'+str(npop)+':'+str(ct)+'.csv'  
            #print 'Running New Ecosim..'
 #           es1=time.time()
 #           runecosim(id)
 #           es1e=time.time()
 #           print 'Finished ecosim in '+str(es1e-es1)+' seconds'
            id=''
 #           print ''
    print 'Loop '+str(ct)+' DONE!'
    print ''
    ct=ct+1
'''            
end=time.time()
print ''
print ''

print 'Finished '+str(totalrun)+' runs'       
print "in " + str(end-start) + " seconds"'''







