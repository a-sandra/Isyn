#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: cdivision=True


import numpy
cimport numpy
cimport cython
cdef class particles:
    cdef public numpy.ndarray particleArray
    cdef public int particleNumber
    cdef public dict coordinateNames
    cdef public double particleCharge
    def __init__(particles self, numpy.ndarray[numpy.float64_t,ndim=2,mode="c"] coordinates,double particleCharge=1.0, coordinateNames={"x":0,"y":1}):
        self.particleArray=coordinates.copy()
        self.particleCharge=particleCharge
        self.particleNumber=coordinates.shape[0]
        self.coordinateNames=coordinateNames
cdef class Cgrid1D:
    cdef public numpy.float64_t left
    cdef public numpy.float64_t right
    cdef public unsigned int nBins
    cdef public numpy.float64_t gridSpacing
    cdef public numpy.ndarray density
    cdef public numpy.ndarray bins
    def __init__(self,double left=0,double right=1.,int nBins=100):
        '''
        expects:
            left=0, float giving left grid border
            right=0, float giving right grid border
            nBins=100, int giving number of grid cells
        '''
        self.left=left
        self.right=right
        self.gridSpacing=1.*(right-left)/nBins
        self.nBins=nBins
        self.density = numpy.zeros(nBins,dtype=numpy.float64)
        self.bins = numpy.arange(left,right,self.gridSpacing)
    cpdef object particles2Grid(Cgrid1D self, particles particlesIn,char[] direction="x"):
        ''' Convenience method
        allows to calculate the densities of particles on a grid.
        expects a string for the direction (direction="x", "y" or something given in the dictionary of 
        '''
        self.points2Grid(particlesIn.particleArray,particlesIn.coordinateNames[direction],pointMass=particlesIn.particleCharge)
    cpdef object points2Grid(Cgrid1D self,numpy.ndarray[numpy.float64_t,ndim=2,mode='c'] particlesIn,direction=0,pointMass=1.):
        '''
        expects:
            a two-dimensional numpy array (in c-mode) of the form 
            [ [p1.c1, p1.c2,..., p1.cn], [p1.c2,p1.c2,..., p2.cn]... [pm.c1,pm.c2,..., pm.cn]]
            If your array is not in C-Mode (for example if you transposed it shortly before) you will likely get an error. In that case you can just
            pass a copy of your original array.
        '''
        cdef double [:,:] tst=particlesIn
        self._points2Grid(&particlesIn[0,0],direction,particlesIn.shape[0],particlesIn.shape[1],pointMass)
    
    cdef void _points2Grid(Cgrid1D self,double* particles, unsigned int direction, unsigned int numberOfParticles, unsigned int dimensions,double pointMass):
        '''
        expects:
            particles: array of coordinates (1-dimensional)
            grid: grid object
            particle coordinates should not exceed grid area
        '''
        cdef numpy.ndarray[numpy.float64_t,ndim=1] tmpgrid=numpy.zeros(self.nBins,dtype=numpy.float64)
    #    cdef double [:] particles=particlesIn
        cdef unsigned int leftIndex =0
        cdef numpy.float64_t binPosition= 0.
        cdef numpy.float64_t invBinWidth=pointMass/self.gridSpacing
        cdef numpy.float64_t left=self.left
        cdef numpy.float64_t right=self.right
        cdef numpy.float64_t gridSpacing=self.gridSpacing
        cdef unsigned int nBins=self.nBins
        cdef unsigned int i = direction

        #for i in range(2*numberOfParticles,):
        # divmod is a python function, so we are faster without it.
        cdef unsigned int totalElements=(numberOfParticles*dimensions)
        while i<totalElements:
            leftIndex=<unsigned int>((particles[i]-left)//gridSpacing) #this is the div part of divmod
            binPosition=((particles[i]-left)%gridSpacing)*invBinWidth#the mod part of divmod
            tmpgrid[(leftIndex % nBins)]+=pointMass-binPosition #interpolation on the neighbouring bins
            tmpgrid[((leftIndex + 1) % nBins)]+=binPosition
            i+=dimensions
        self.density=tmpgrid
    
cdef class CGrid2D:
    cdef public numpy.float64_t xmin,xmax,ymin,ymax
    cdef public unsigned int nBinsX,nBinsY
    cdef public numpy.float64_t gridSpacingX,gridSpacingY
    cdef public numpy.ndarray density
    cdef public numpy.ndarray bins
    def __init__(self,double xmin=0,double xmax=1.,double ymin=0, double ymax=1,int nBinsX=100,nBinsY=100):
        '''
        expects:
            left=0, float giving left grid border
            right=0, float giving right grid border
            nBins=100, int giving number of grid cells
        '''
        self.xmin=xmin
        self.xmax=xmax
        self.ymin=ymin
        self.ymax=ymax
        (self.gridSpacingX,self.gridSpacingY)=(1.*(xmax-xmin)/nBinsX,1.*(ymax-ymin)/nBinsY)
        self.nBinsX=nBinsX
        self.nBinsY=nBinsY
        self.density = numpy.zeros((nBinsX,nBinsY),dtype=numpy.float64)
       # self.bins = numpy.arange((),self.gridSpacing)
    def points2Grid(self,numpy.ndarray[numpy.float64_t,ndim=2,mode='c'] particlesIn,double pointMass=1.,unsigned int dimensions=2):
        cdef tmpgrid= numpy.zeros((self.nBinsX,self.nBinsY),dtype=numpy.float64)
        cdef numpy.float64_t invBinWidthX=pointMass**0.5/self.gridSpacingX
        cdef numpy.float64_t invBinWidthY=pointMass**0.5/self.gridSpacingY
        for particle in particlesIn:
            particleX=particle[0]
            particleY=particle[1]
            leftIndexX=<unsigned int>((particleX-self.xmin)//self.gridSpacingX)
            leftIndexY=<unsigned int>((particleY-self.ymin)//self.gridSpacingY)
            binPositionX=((particleX-self.xmin)%self.gridSpacingX)*invBinWidthX
            binPositionY=((particleY-self.ymin)%self.gridSpacingY)*invBinWidthY
            tmpgrid[leftIndexX%self.nBinsX][leftIndexY%self.nBinsY]+=(pointMass**0.5-binPositionX)*(pointMass**0.5-binPositionY)
            tmpgrid[(leftIndexX+1)%self.nBinsX][leftIndexY%self.nBinsY]+=(binPositionX)*(pointMass**0.5-binPositionY)
            tmpgrid[(leftIndexX+1)%self.nBinsX][(leftIndexY+1)%self.nBinsY]+=(binPositionX)*(binPositionY)
            tmpgrid[leftIndexX%self.nBinsX][(leftIndexY+1)%self.nBinsY]+=(pointMass**0.5-binPositionX)*(binPositionY)
        print tmpgrid
        self.density=tmpgrid