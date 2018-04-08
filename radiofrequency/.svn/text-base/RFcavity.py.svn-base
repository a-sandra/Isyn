import numpy as np
import math
from physicalconstant.PhysicalConstant import *
from beamofparticle.Particle import Particle
from accelerator.Accelerator import Accelerator

class RFcavity(object):
    """ Class RFcavity defined by:
    the harmonic number
    the RF voltage
    """
#    cdef public double harmonic
#    cdef public double voltage
#    cdef public double stablephase
#    cdef public double eta
#    cdef object accObject
#    cdef double gammaParticle
    
    def __init__(self, accObject):
        self.harmonic = 1
        self.voltage = 1
        self.stablephase = 0 # Would be nice to have it automatically
        self.eta = 0.001 # Temporary
        self.accObject = accObject
    
    def setHarmonic(self, valueh):
        self.harmonic = valueh
    
    def getHarmonic(self):
        return self.harmonic
    
    def setRFvoltage(self, valuev):
        self.voltage = valuev

    def getRFvoltage(self):
        return self.voltage
    
    def setStablePhase(self, valuephis):
        self.stablephase = valuephis
        
    def getStablePhase(self):
        return self.stablephase
    
    def setEta(self,valueEta):
        self.eta = valueEta
    
    def getEta(self):
        return self.eta

    def rfkick(self, beamObject):
#        print self.eta, beamObject.qCharge, self.voltage, beamObject.beta, self.stablephase, self.harmonic
        beamObject.particleArray[:,1] += beamObject.qCharge *e* self.voltage/(beamObject.beta**2*beamObject.TotalEnergy)*(np.sin(beamObject.particleArray[:,0]) - np.sin(self.stablephase))
        eta=lambda dpp: (1/self.accObject.gammat)**2-(1/(beamObject.gamma*(dpp*beamObject.beta**2 + 1)))**2
        # def eta(dpp): return (1/self.accObject.gammat)**2-(1/(beamObject.gamma*(dpp*beamObject.beta**2 + 1)))**2
        beamObject.particleArray[:,0] += 2.0*math.pi * self.harmonic* eta(beamObject.particleArray[:,1]) * beamObject.particleArray[:,1]
        
    def track(self, beamObject):
        self.updateStablePhase(beamObject)
        self.rfkick(beamObject)
    
    def updateStablePhase(self,beamObject):
        if self.accObject.eta*math.cos(self.stablephase) >= 0:
            print 'Transition'
            self.phasetemp = self.stablephase
            self.stablephase = math.pi - self.stablephase
            beamObject.particleArray[:,0]=beamObject.particleArray[:,0]+(math.pi-2*self.phasetemp)

    def printRFparameters(self):
        print "harmonic number =", self.harmonic
        print "RF voltage [V]= ", self.voltage
        print "Stable phase in [rad] = ", self.stablephase