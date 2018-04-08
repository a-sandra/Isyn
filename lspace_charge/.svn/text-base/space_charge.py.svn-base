#import os
#import sys
from physicalconstant.PhysicalConstant import *
from accelerator.Accelerator import Accelerator
from beamofparticle.Beam import Beam
from radiofrequency.RFcavity import RFcavity
import numpy as np
from part2grid.particles2grid import *
import matplotlib.pyplot as plt

class spaceChargeNode(object):
    """ Class for longitudinal space charge kick"""
    
    def __init__(self, accObject):
        self.macroparticle = 1000000 # to be change to get from Nparticle acce. object
        self.nbins = 128
        self.ratio = 5
        self.eta = 0.001
        self.totalnumberOfParticle = 2.0e13
        self.offsetphi = 0.0
        self.left = 0.0
        self.right = 0.0
        self.gridSpacing = 1.*(self.right-self.left)/self.nbins
        self.harmonic = 5 # Here, we need to give the harmonic to the space charge module. Tricky. Temporary solution
        self.Tsamp = 0.001
        self.charge = 1.0* e
        self.accObject = accObject
        self.circum = 1.
        self.g0 = 3.84
        self.initialDeltaPhi = 0.48 # has to be changed
    
    def setMacroparticle(self, nbMacro): # not needed
        self.macroparticle = nbMacro
    
    def setNumberParticleinBunch(self, nbpart):
        self.totalnumberOfParticle = nbpart
    
    def setNumberBins(self, bins):
        self.nbins = bins
    
    def getNumberBins(self):
        return self.nbins
    
    def setRatio(self, rat):
        self.ratio = rat
    
    def setEta(self,valueEta):
        self.eta = valueEta
        
    def track(self, beamObject):
        self.space_charge_kick(beamObject)
        
    def ParticleShape(self,x):
        return max(1-abs(x),0)
    
    def spaceChargeImpedance(self, beamObject, w):
        return (2*math.pi)*w*mu0*self.circum/(4*math.pi*(beamObject.beta*beamObject.gamma)**2)*self.g0
    
    def phi2z(self,har):
        self.harmonic = har
    
    def updateStablePhase(self,beamObject):
        pass
        
    def ScVoltage(self, beamObject, x):
        return 3.0/4.0*(self.charge*self.g0*self.totalnumberOfParticle*self.harmonic**2)/(eps0*beamObject.gamma**2*self.initialDeltaPhi**3*self.accObject.radius)*(-self.harmonic/self.accObject.radius*x)

# --------- Main routine for space charge kick ----

    def space_charge_kick(self, beamObject):
        self.circum = 2*math.pi*self.accObject.radius
        omegar = beamObject.beta * c * self.accObject.radius
        self.charge = beamObject.qCharge * e
        chargePerParticle = self.totalnumberOfParticle/beamObject.nParticle
        self.g0 = 0.5 +2.0*math.log(self.ratio)
        self.left = -(math.pi - self.offsetphi)*self.accObject.radius/self.harmonic
        self.right= (math.pi - self.offsetphi)*self.accObject.radius/self.harmonic
        self.gridSpacing=1.*(self.right-self.left)/self.nbins
        ipart = beamObject.beta * c * self.charge * self.totalnumberOfParticle/beamObject.nParticle # in m/s x C/m current in a grid cell
        coords=np.array([ (beamObject.particleArray[:,0]-self.offsetphi)*self.accObject.radius/self.harmonic,np.zeros(beamObject.nParticle)]).T
        grid = Cgrid1D(self.left,self.right,self.nbins)
        grid.left = self.left
        grid.right = self.right
        grid.gridSpacing = self.gridSpacing
        grid.nBins = self.nbins
        grid.points2Grid(coords.copy(),direction = 0 , pointMass = ipart)
        #print (np.sum(grid.density))
        self.Tsamp = self.gridSpacing/(beamObject.beta * c)
        freq = np.fft.fftfreq(self.nbins,d=self.Tsamp)
        fftIntensity = np.fft.fft(grid.density,self.nbins)
        productZsc2Intensity = -1J*self.spaceChargeImpedance(beamObject,freq)*fftIntensity
        fftInverse = np.fft.ifft(productZsc2Intensity).real
        grid.density = fftInverse.copy()
        coordinates=coords.copy()
        grid.grid2Points(coordinates, directionIn = 0 , pointMass = beamObject.qCharge)
        # coordinates[:,1] is the space charge voltage
        beamObject.particleArray[:,1] += 1/(beamObject.beta**2*beamObject.TotalEnergy) * beamObject.qCharge *e* coordinates[:,1]/(2*math.pi*self.accObject.radius)*beamObject.beta * c
        eta=lambda dpp: (1/self.accObject.gammat)**2-(1/(beamObject.gamma*(dpp*beamObject.beta**2 + 1)))**2
        beamObject.particleArray[:,0] += 2.0*math.pi * self.harmonic* eta(beamObject.particleArray[:,1]) * beamObject.particleArray[:,1]