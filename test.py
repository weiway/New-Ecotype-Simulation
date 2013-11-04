import time
import os

def runecosim(runid,version,path_to_fasta,path_to_tree):
    path='ecosim'+version
    outname='out_'+version+'_'+str(runid)
    os.chdir(path)
    cwd=os.getcwd()

    start=time.time()
    cmd='java -jar bin/EcoSim.jar -n ' + path_to_fasta + path_to_tree + outname
    os.system(cmd)
    end=time.time()

    os.chdir('../')
    print str(runid),"took",str(end-start),"seconds" 
    return end-start

stuff = [('cleaned_1.fas ','cleaned_tree_1.txt '),
        ('cleaned_2.fas ','cleaned_tree_2.txt '),
        ('cleaned_3.fas ','cleaned_tree_3.txt '),
        ('cleaned_4.fas ','cleaned_tree_4.txt '),
        ('cleaned_5.fas ','cleaned_tree_5.txt ')]

ite = 1
for i in stuff:
    n=1
    t2l=[]
    while n<21:
        print 'running round '+str(n)
        t2=runecosim(n,'2.14','../sequences/'+i[0],'../sequences/'+i[1])
        t2l.append(t2)
        n+=1


    print
    print '2.14 time usage:'
    for t in t2l:
        print str(t) + ' secs'

    avg2 = sum(t2l)/float(len(t2l))

    print 
    print '2.14 avg time: ' + str(avg2) + ' secs'
    ite+=1 







