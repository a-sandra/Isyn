import math

class Accelerator_Info(object):
    """ Summaries the info of the accelerator and beam"""
    def __init__(self,accObject,beamObject):
        self.eta = math.pow(1/(accObject.gammat), 2) - math.pow(1/(beamObject.gamma),2) #See if we cannot do in the other way
        print "gamma = ", beamObject.gamma
        print "gamma transition = ", accObject.gammat
        print "eta = ", self.eta
        print "Bdot = ", accObject.bdot