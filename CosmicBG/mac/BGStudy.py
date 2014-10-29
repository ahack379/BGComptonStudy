import sys

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
gSystem.Load("libBGComptonStudy_CosmicBG")
from ROOT import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()

#for x in xrange(len(sys.argv)-6):
#    my_proc.add_input_file(sys.argv[x+1])

# Set input root file
my_proc.add_input_file(sys.argv[1])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify input TDirectory name if given
if len(sys.argv) > 2:

    my_proc.set_input_rootdir(sys.argv[2])

# Specify output root file name
my_proc.set_ana_output_file("out.root");

# Attach a template process
my_proc.add_process(fmwk.BGShowerInfo());

print
print  "Finished configuring ana_processor. Start event loop!"
print

# Let's run it.
my_proc.run();

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

