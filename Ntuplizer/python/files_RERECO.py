redirector='root://gfe02.grid.hep.ph.ic.ac.uk:1097/'

BToKee_Seed2p0 = []
base = redirector+'/store/user/bainbrid/BToKee_Pythia/LowPtEle_2GeV/180508_130015/'
for index1 in range(7) :
    for index2 in range(1000) : 
        index = index1*1000 + index2
        if index == 0 or index > 6947 : continue 
        BToKee_Seed2p0.append( base + '000{:.0f}/'.format(index1) + 'RECOSIM_{:.0f}.root'.format(index) )

BToKee_Seed1p0 = []
base = redirector+'/store/user/bainbrid/BToKee_Pythia/LowPtEle_1GeV/180508_150737/'
for index1 in range(7) :
    for index2 in range(1000) : 
        index = index1*1000 + index2
        if index == 0 or index > 6947 : continue 
        BToKee_Seed1p0.append( base + '000{:.0f}/'.format(index1) + 'RECOSIM_{:.0f}.root'.format(index) )

BToKee_Seed0p5 = []
base = redirector+'/store/user/bainbrid/BToKee_Pythia/LowPtEle_0p5GeV/180514_204455/'
for index1 in range(7) : 
    for index2 in range(1000) : 
        if index2 in [6014] : continue # failed jobs
        index = index1*1000 + index2
        if index == 0 or index > 6947 : continue 
        BToKee_Seed0p5.append( base + '000{:.0f}/'.format(index1) + 'RECOSIM_{:.0f}.root'.format(index) )

BToKee_Gsf0p5 = []
base = redirector+'/store/user/bainbrid/BToKee_Pythia/EleReco-Gsf0p5/180606_221056/'
for index1 in range(7) : 
    for index2 in range(1000) : 
        if index2 in [] : continue # failed jobs (none)
        index = index1*1000 + index2
        if index == 0 or index > 6947 : continue 
        BToKee_Gsf0p5.append( base + '000{:.0f}/'.format(index1) + 'EleReco-Gsf0p5_{:.0f}.root'.format(index) )

BToKmm_Seed2p0 = []
base = redirector+'/store/user/bainbrid/BToKmm_Pythia/LowPtMuon_2GeV/180523_132754/'
for index1 in range(9) :
    for index2 in range(1000) : 
        if index2 in [8855,] : continue # failed jobs #@@ if problems, add 7,8,9,10 ? 
        index = index1*1000 + index2
        if index == 0 or index > 9525 : continue 
        BToKmm_Seed2p0.append( base + '000{:.0f}/'.format(index1) + 'RECOSIM_{:.0f}.root'.format(index) )

BToKmm_Seed0p5 = []
base = redirector+'/store/user/bainbrid/BToKmm_Pythia/LowPtMuon_0p5GeV/180514_213021/'
for index1 in range(9) :
    for index2 in range(1000) : 
        if index2 in [1,2,3,4,5,21,82,157,238,319,396,463,538,] : continue # failed jobs
        index = index1*1000 + index2
        if index == 0 or index > 9525 : continue 
        BToKmm_Seed0p5.append( base + '000{:.0f}/'.format(index1) + 'RECOSIM_{:.0f}.root'.format(index) )

if __name__ == '__main__' :
    # B->Kee
    for afile in BToKee_Seed1p0 : print afile; print 
    for afile in BToKee_Seed2p0 : print afile; print 
    for afile in BToKee_Seed0p5 : print afile; print 
    for afile in BToKee_Gsf0p5  : print afile; print 
    # B->Kmm
    for afile in BToKmm_Seed2p0 : print afile; print 
    for afile in BToKmm_Seed0p5 : print afile; print 
