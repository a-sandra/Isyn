import numpy as np
import math
import random
from Particle import Particle
import matplotlib.pyplot as plt
from physicalconstant.PhysicalConstant import *

class Beam(object):
    """ Class Beam defined by:
    - Its relativistic gamma
    - Its rest energy
    - Kind of distribution
    - Number of particles """
    
    def __init__(self): # Method constructor
        self.gamma = 1.0
        self.Erest = 938272046.0
        self.qCharge = 1.0
        self.nParticle = 100000
        self.bunchOfParticle = []
        self.particleArray=np.zeros((self.nParticle,2),dtype="float64")
        self.beta = math.sqrt(1- 1/math.pow(self.gamma, 2))
        self.distributionKind = "Parabolic"
        self.givendistribution = []
        
    def setGamma(self,valuegamma):
        self.gamma = valuegamma
        self.updateBeamParameters()
    
    def getGamma(self):
        return self.gamma
    
    def setRestEnergy(self, valueE0):
        self.Erest = valueE0
    
    def getRestEnergy(self):
        return self.Erest
    
    def setNParticle(self, valueNb):
        self.nParticle = valueNb
        self.particleArray=np.zeros((self.nParticle,2),dtype="float64")
    
    def getNParticle(self):
        return self.nParticle
    
    def setCharge(self, valuecharge):
        self.qCharge = valuecharge
        
    def updateBeamParameters(self):
        self.beta = math.sqrt(1- 1/math.pow(self.gamma, 2))
        self.TotalEnergy = self.gamma*self.Erest*e
        
    def distributionFromFile(self, filename):
        self.givendistribution = np.loadtxt(filename)
        for i in range(self.nParticle):
            self.bunchOfParticle.append(Particle(self.givendistribution[i][0], self.givendistribution[i][1]))
            print self.givendistribution[i][0]
        self.particleArray=self.givendistribution

    def parabolicDistribution(self, initialDeltaPhi, initialdp, offsetphi = 0):
        self.distributionKind = "Parabolic"
        random.seed(3000)
        for i in range(self.nParticle):
            u = random.uniform(0,1)
            v = random.uniform(0,1)
            phi = offsetphi + initialDeltaPhi * math.sqrt(1.0 - math.pow(1.0 - u, 2.0/3.0))*math.cos(2*math.pi*v)
            dpi = initialdp * math.sin(2.0*math.pi*v)*math.sqrt(1.0 - math.pow(1.0- u, 2.0/3.0))
            self.bunchOfParticle.append(Particle(phi, dpi))
            self.particleArray[i,0]=phi
            self.particleArray[i,1]=dpi

    def getParticles(self):  
        return self.particleArray.T
    
    def plotBunch(self):
        plt.plot(self.particleArray[:,0], self.particleArray[:,1],"r.")
        plt.show()
    
    def plotLongitudinalDistribution(self):
        phi,dp = self.getParticles()
        plt.hist(phi,bins = 100, color = 'green' )
        plt.show()
    
    def plotEnergyDistribution(self):
        phi,dp = self.getParticles()
        plt.hist(dp,bins = 100, color = 'green' )
        plt.show()

    def getRmsValues(self):
        phi, dp = self.getParticles()
        phi_rms = np.std(phi)
        dp_rms = np.std(dp)
        return phi_rms , dp_rms

    def getMean(self):
        phi, dp = self.getParticles()
        phi_mean = np.mean(phi)
        dp_mean = np.mean(dp)
        return phi_mean , dp_mean

    def printBeamParameters(self):
        """ Print to screen a summary of the attributes of the object Beam """
        print "beta=", self.beta
        print "gamma=", self.gamma
        print "Rest Energy=", self.Erest
        print "Total Energy = ", self.TotalEnergy
        print "Charge of particle = ", self.qCharge
        print "Number of particle = ", self.nParticle
        print "Type of distribution =", self.distributionKind
