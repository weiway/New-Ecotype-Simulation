










def runecosim(runid,version,path_to_fasta,path_to_tree):

    import time
    import os
    path='./ecosim'+version
    outname='out_'+version+'_'+str(runid)
    os.chdir(path)
    cwd=os.getcwd()

    start=time.time()
    cmd='java -jar bin/EcoSim.jar -n ' + path_to_fasta+ path_to_tree + outname
    os.system(cmd)
    end=time.time()

    os.chdir('../')
    return end-start

n=1
t1l=[]
t2l=[]
while n<21:
    print 'running round '+str(n)
    t1=runecosim(n,'2.13c','../sequences/AP_Edited.fas',' ../sequences/APtree_Edited.nwk ')
    t2=runecosim(n,'2.14','../sequences/AP_Edited.fas',' ../sequences/APtree_Edited.nwk ')
    t1l.append(t1)
    t2l.append(t2)
    n+=1

print '' 
print '2.13c time usage:'
for t in t1l:
    print str(t) + ' secs'
avg1 = sum(t1l)/len(t1l)


print ''
print '2.14 time usage:'
for t in t2l:
    print str(t) + ' secs'

avg2 = sum(t1l)/len(t1l)

print ''
print '2.13c avg time: ' + str(avg1) + ' secs'
print '2.14 avg time: ' + str(avg2) + ' secs'







