from __future__ import division

cdef class Particle:
    cpdef public double phia,dpa
    def __init__(self,phia,dpa):
        self.phia=phia
        self.dpa=dpa
        pass