import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

import os, ROOT
from ROOT import gSystem
gSystem.Load("libLArLite_Base")
gSystem.Load("libLArLite_Analysis")
gSystem.Load("libLArLite_LArUtil")
gSystem.Load("libBGComptonStudy_CosmicBG")
gSystem.Load("libMCPartInfo_MCPartGetter")
gSystem.Load("libBasicTool_GeoAlgo")

from ROOT import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

for x in xrange(len(sys.argv)-1):
    my_proc.add_input_file(sys.argv[x+1])

# Set input root file
#my_proc.add_input_file(sys.argv[1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify input TDirectory name if given
if len(sys.argv) > 2:

    my_proc.set_input_rootdir(sys.argv[2])

# Specify output root file name
my_proc.set_ana_output_file("out.root");

# Attach a template process

bg = fmwk.CosmicsBackground();
bg.setECut(0.01)#GeV

# Vector of PDGs to search
pdgs    = ROOT.vector('int')()
pdgs.push_back(11)
pdgs.push_back(13)
pdgs.push_back(-13)

mcgetter = fmwk.MCgetter()
mcgetter.getAllPDGs(pdgs)
# Energy cut: If PDG match && E > _ECut [GeV] then add particle
mcgetter.SetECut(0.01)
mcgetter.SetVerbose(False)

# Tell module what PDGs to search for
bg.SetMCgetter(mcgetter)


my_proc.add_process(bg);

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run(0,18);

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

