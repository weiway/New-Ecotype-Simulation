
def produce_sequences(omega, sigma, npop, drift, nu, lengthseq):
    #The produce sequences input file
    acinas = open("./producesequences/acinas.dat",'w')
    acinas.write(str(omega) + '      omega\n')
    acinas.write(str(sigma) + '      sigma\n')
    acinas.write(str(npop) + '        npop\n')
    acinas.write(drift+ '   drift\n')
    acinas.write(nu + '       nu\n')
    import random
    acinas.write(str(random.randint(0,2**31-1)) + '           iii (random number seed)\n')
    acinas.write(lengthseq + '      lengthseq (after deleting gaps, etc.)')
    acinas.close()
    import os
    import platform
    
    #Run the produce sequences script from producesequences
    os.chdir('./producesequences/')
    if platform.system() == 'Windows':
        os.system('acinaspolz.exe')
    else:
        import shutil
        import os.path

        shutil.copyfile('acinas.dat','./acinaspolz.app/Contents/Resources/drive_c/Program Files/acinaspolz/acinas.dat')

        os.system('open acinaspolz.app')
        os.chdir('./acinaspolz.app/Contents/Resources/drive_c/Program Files/acinaspolz/')
        
        while os.path.isfile('produceseqs.dat')==False:
            import time
            time.sleep(1)
            print 'generating... please wait...'
        #Aaron: CHANGE PATH HERE
        shutil.copyfile('produceseqs.dat','/Users/Wei/Desktop/produceseq/producesequences/produceseqs.dat')
        os.remove('produceseqs.dat')
    os.chdir('../')

