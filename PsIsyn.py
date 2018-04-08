import os
import sys
sys.path.append("~/Documents/workspace/Isyn/")
from accelerator.Accelerator import Accelerator
from accelerator.Accelerator_Info import Accelerator_Info
from beamofparticle.Beam import Beam
from radiofrequency.RFcavity import RFcavity
from lspace_charge.space_charge import spaceChargeNode
import cProfile
# Version of Matplotlib
#print "Version of matplotlib", mpl.__version__
"""
os.system('rm output_file/*')
os.system('mkdir output_file')"""

#--------------- Declare the Object Accelerator ------------#

ps = Accelerator()
ps.setRadius(5000.0)
ps.gammat = 6.1 #ps.bdot = 4.0
ps.bdot = 0.0
ps.gammadot = 0.0

ps.setParameterFilename('ps')
ps.dumpDistributionEveryNturn(100000)

# --------------- Declare the Object Beam ------------#

beam = Beam()
beam.setGamma(2.12)
beam.setRestEnergy(938272046.0)
beam.setCharge(1.0)
beam.setNParticle(1000000) #beam.distributionFromFile('AccInitDistribution.dat') #beam.parabolicDistribution(0.396017, 0.00161535, 0.943477)
beam.parabolicDistribution(0.48, 0.001418, 0.0)
#beam.parabolicDistribution(0.509, 0.001418, 0.0)
#beam.plotBunch()

# --------------- Get the Object Accelerator_Info ------------#

info = Accelerator_Info(ps, beam)

# --------------- Declare the Object Rfcavity ------------#

rf = RFcavity(ps)
rf.setHarmonic(5)    #rf.setRFvoltage(280000.0)
rf.setRFvoltage(60000.0)
rf.setStablePhase(0.0)    #rf.setStablePhase(0.943477)
ps.addElement(rf)

# --------------- longitudinal space charge node ------------#

sc_node = spaceChargeNode(ps)
ps.addElement(sc_node)
sc_node.setMacroparticle(1000000) # redondant, not necessary, i will have to find something else.
sc_node.setNumberBins(128)
sc_node.setRatio(5)
sc_node.phi2z(rf.getHarmonic()) # This is temporary, because next month, the code will be slightly reshuffled by
#introducing the grid with the RFcavity.

# --------------- Tracking ------------#

cProfile.run('ps.track(beam, 1)','restats')
import pstats
p = pstats.Stats('restats')
p.sort_stats('cumulative').print_stats(10)

#beam.plotBunch()

# ------- Additional stuff -----#
"""
os.system('mv ps* output_file')
os.system('mv *.png output_file/')
"""
