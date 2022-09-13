#!/usr/bin/env python

# rspecies - creates rspecies file from rreacs file
# checkspec - checks to make sure that there is a production and destruction reaction for every species

######################################################################
# Create the rspecies file to correspond to the rreacs file
def rspecies(rreacs):
    print( "rspecies")
    fin = open(rreacs,'r')

    species = []

    for line in fin:
        if line[0] == '#': continue

        mols = []
        words = line[:86].split()

        for spec in words:
            if spec in species: continue
            try:
                float(spec)
                continue
            except ValueError:
                species.append(spec)

    fin.close()

    notmols = ['CRP', 'PHOTON', 'XRAY', 'TEMP', 'C-RAY', 'LYAPHOTON']
    rspecies = rreacs.replace('reacs','species')
    fout = open(rspecies, 'w')    
    fout.write("# %d\n" % (len(species) - len(notmols)))


    for mol in species:
        if mol in notmols: continue
        fout.write("%s\n" % mol)

    for mol in notmols:
        fout.write("# %s\n" % mol)

    fout.close()
    
######################################################################
# for every species in rspecies, check to make sure that there is a production and
# destruction reaction in rreacs
def checkspec(rreacs):
    print( "checkspec")
    rspecies = rreacs.replace('reacs','species')

    species = {}
    fspecies = open(rspecies, 'r')
    fspecies.readline()
    for line in fspecies:
        if line[0] == '#':
            species[line.split()[1]] = ''
        else:
            species[line.split()[0]] = ''
    fspecies.close()
    
    freacs = open(rreacs, 'r')
    for line in freacs:
        if line[0] == '#': continue
        reactants = line[6:32].split()
        products = line[32:86].split()

        for val in reactants: species[val] += 'r'
        for val in products: species[val] += 'p'

    freacs.close()

    for spec in species.keys():
        if species[spec].find('r') == -1 or species[spec].find('p') == -1:
            print( spec)

######################################################################
# for every neutral species in rspecies, check to make sure that there are
# adsorption and desorption reactions
def checkgr(rreacs):
    print( "checkgr")
    rspecies = rreacs.replace('reacs','species')

    species = {}
    fspecies = open(rspecies, 'r')
    for line in fspecies:
        if line[0] == '#': 
            continue
        spec = line.split()[0]
        if '-' not in spec and '+' not in spec and '(gr)' not in spec:
            species[line.split()[0]] = False
    fspecies.close()
    
    freacs = open(rreacs, 'r')
    for line in freacs:
        if line[0] == '#': continue
        reactants = line[6:32].split()

        if 'GRAIN-' in line: continue
        elif 'GRAIN+' in line: continue
        elif 'GRAIN0' in line: continue
        elif 'GRAIN' in line: 
            for val in reactants:
                if val != 'GRAIN': species[val] = True

    freacs.close()

    for spec in species.keys():
        if species[spec] == False:
            print( spec)

######################################################################

def main():
    rreacs = 'rreacs_osu2009_umfmt.dat'
    print( "call rspecies")
    rspecies(rreacs)    # create rspecies file
    checkspec(rreacs)   # check if all species have formation/destruction reaction 
    #checkgr(rreacs)     # check if all neutral species have grain reactions
    
main()
