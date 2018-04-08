# ----- Import the useful Python packages ----- #
import numpy as np
import math
import matplotlib.pyplot as plt
#from sympy.solvers import solve

# ----- Import the customed packages ----- #
from physicalconstant.PhysicalConstant import *
import misc.Header

# ----- Class Accelerator----- #

class Accelerator(object):
    
    def __init__(self):
        self.elements = []
        self.radius = 172.4603
        self.gammat = 8.9
        self.gammadot = 65.6221
        self.bdot = 4.0
        self.eta = 0.001
        self.time = 0.0
        self.phis=0.0
        self.parameterfilename = 'sis100Default'
        self.nturndump = 1000
        self.i =1
    
    def setParameterFilename(self,filename):
        self.parameterfilename = filename + ".dat"
        self.parameterfile = file(self.parameterfilename,'w')
        self.parameterfile.write('Time (s)  TotalEnergy (SI)  gamma  beta  eta  Revolution time (s) phi_s \n')
        self.beamfilenameRms = self.parameterfilename + '_Rms.dat'
        self.beamfileRms = file(self.beamfilenameRms,'w')
        
    def WriteToFileBeamParameter(self,beamObject):
        self.parameterfile.write('%.10f %.10g %.10f %.10f %.10f %.10f %.5f \n' %(self.time, beamObject.TotalEnergy, beamObject.gamma, beamObject.beta, self.eta,2*math.pi*self.radius/(beamObject.beta*c), self.phis))
        phi_rms , dp_rms = beamObject.getRmsValues()
        phi_mean , dp_mean = beamObject.getMean()
        self.beamfileRms.write('%.6f %.6f %.6f %.6f \n ' %(phi_rms , dp_rms, phi_mean , dp_mean))
        
    def dumpDistributionEveryNturn(self, nbturn):
        self.nturndump = nbturn

    def setBdot(self,valueBdot):
        self.bdot = valueBdot
    
    def gammadot(self,valuegammadot):
        self.gammadot = valuegammadot

    def addElement(self, element):
        self.elements.append(element)
        print self.elements
    
    def setRadius(self,valueRadius):
        self.radius = valueRadius
    
    def getRadius(self):
        return self.radius
    
    def setGammaTr(self,valueGammaTr):
        self.gammat = valueGammaTr
    
    def computeEta(self,beamObject):
        self.eta = math.pow(1/(self.gammat), 2) - math.pow(1/(beamObject.gamma),2)
    
    """ Ideally updateParameter and updatephis are not relevant for a space charge node"""
    """ updateparameter could be in rfcsvity as well as updatephis"""
    def track(self, beamObject, numberOfTurn):
        self.computeEta(beamObject)
        for i in np.arange(numberOfTurn):
            for element in self.elements:
                element.setEta(self.eta)
                element.track(beamObject)
                self.updateParameter(beamObject) #Later with impedance etc.. this command should be element.updateParam.. in which for a RF cavity the phase is updated
                element.updateStablePhase(beamObject)
#                self.phis = element.getStablePhase() # to be change, does not work with SC kick
            self.WriteToFileBeamParameter(beamObject)
#            self.diagnosticBeam(i, beamObject)
#               .element.update() # to implement later

    def updateParameter(self,beamObject):
        self.time = self.time + 2*math.pi*self.radius/(beamObject.beta*c)
        DeltaSynchrotronEnergy = self.gammadot*beamObject.Erest*e*2*math.pi*self.radius/(beamObject.beta*c) # in eV
        beamObject.TotalEnergy = beamObject.TotalEnergy + DeltaSynchrotronEnergy # in eV
        beamObject.gamma = beamObject.TotalEnergy/(beamObject.Erest*e) #in eV
        beamObject.beta = math.sqrt( 1 - 1/math.pow(beamObject.gamma, 2))
        self.computeEta(beamObject)
        

    def diagnosticBeam(self,nbofturn, beamObject):
        if (math.fmod(nbofturn, self.nturndump)==0):
            self.filenamedistr = open(self.parameterfilename + 'DumpDistr_'+str(nbofturn)+'.dat','w')
            phi,dp = beamObject.getParticles()
            np.savetxt(self.filenamedistr, np.array([phi,dp]).T)
            self.filenamedistr.close()
            plt.plot(phi,dp,"r.")
            plt.xlabel('phase')
            plt.ylabel('momentum spread')
            plt.axis([0, math.pi, -0.005, 0.005])
            plt.title(str(nbofturn)+'turns - eta = '+str(self.eta))
            plt.savefig('fig%05d.png' %self.i)
            plt.close()
            self.i=self.i+1

