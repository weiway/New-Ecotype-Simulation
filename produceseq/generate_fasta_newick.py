#This is working somewhat correctly
#from adaptml.inputconverter import fred_inputconverter

#from runrealecotypes import *

import random

'''def generate(input_file, iso_likelihood, path_to_PhyML_exe):
    try:
        input = open(input_file,'r')
    except IOError:
        print 'Cannot open ' + input_file
    fasta = open('./sequences/fasta.fas','w')
    at_actual_data = 0
    for line in input:
        if at_actual_data == 1:
            fasta.write('\n' + line.strip())
        elif line[0] == '>':
            at_actual_data=1
            fasta.write(line.strip())
    input.close()
    fasta.close()
    fred_inputconverter.convert('./sequences/fasta.fas', iso_likelihood)
    #Now, make the Newick tree.
    import os
    #os.system('cp sequences/amlfredinterleaved.phy output/phyml_current_run_interleaved_2000.fasta')
    os.system(path_to_PhyML_exe +
        ' -i ./sequences/amlfredinterleaved.phy -b 0')
    import shutil
    shutil.copyfile('./sequences/amlfredinterleaved.phy_phyml_tree.txt','./sequences/amlfredinterleaved.phy_phyml_tree_original.txt')
 '''   
def quick_generate(input_file, iso_likelihood, path_to_FastTree):
    try:
        input = open(input_file,'r')
    except IOError:
        print 'Cannot open ' + input_file
    fasta = open('./sequences/fasta.fas','w')
    at_actual_data = 0
    for line in input:
        if at_actual_data == 1:
            fasta.write('\n' + line.strip())
        elif line[0] == '>':
            at_actual_data=1
            fasta.write(line.strip())
    input.close()
    fasta.close()
    #fred_inputconverter.convert('./sequences/fasta.fas', iso_likelihood)
    #Now, make the Newick tree.
    ###
    ###  WE LIE FOR THE NAME OF THE FILE
    ###
    import os
    os.system(path_to_FastTree + ' ' + input_file + ' > ./sequences/fastree.tree.txt')
    import shutil
    shutil.copyfile('./sequences/fastree.tree.txt','./sequences/amlfredinterleaved.phy_phyml_tree_original.txt')
    shutil.copyfile('./sequences/fastree.tree.txt','./sequences/amlfredinterleaved.phy_phyml_tree.txt')
    shutil.copyfile(input_file, './sequences/amlfredfasta.fas')

'''def convert_fasta_and_tree(gamma):
    replaced_strain_names = fred_inputconverter.convert('./sequences/fasta.fas', gamma)
    #return list of new, replaced strain names
    f = open('./sequences/amlfredinterleaved.phy_phyml_tree_original.txt','r')
    treestring = ''
    for line in f.readlines():
        treestring += line
    for strain_name in replaced_strain_names:
        treestring = treestring.replace(getold(strain_name),strain_name)
    f = open('./sequences/amlfredinterleaved.phy_phyml_tree.txt','w')
    f.write(treestring)
    f.close()

def getold(strain_name):
    if strain_name[0] == 'A': return 'B' + strain_name[1:]
    elif strain_name[0] == 'B': return 'A' + strain_name[1:]'''
            
