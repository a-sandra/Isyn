#cython: boundscheck=False
#cython: nonecheck=False
#cython: wraparound=False
#cython: cdivision=True
#cython: embedsignature=True
#True


import numpy
cimport numpy
cimport cython
cdef class particles:
    cdef public numpy.ndarray particleArray
    """Array of particle coordinates """
    cdef public int particleNumber
    """Number of particles"""
    cdef public dict coordinateNames
    """dictionary of coordinate names"""
    cdef public double particleCharge
    """Charge of one particle"""
    cdef public unsigned int dimensions
    """Number of coordinates defining a particle"""
    def __init__(particles self, numpy.ndarray[numpy.float64_t,ndim=2,mode="c"] coordinates,double particleCharge=1.0, coordinateNames={"x":0,"y":1}):
        self.particleArray=coordinates.copy()
        self.particleCharge=particleCharge
        self.particleNumber=coordinates.shape[0]
        self.dimensions=coordinates.shape[1]
        self.coordinateNames=coordinateNames
    cpdef object propagate(self,numpy.float64_t dt,char[] momentumDir="y", char[] locationDir="x"):
        cdef double [:,:] particleArray=self.particleArray
        self._propagate(&particleArray[0,0],dt,momentumDir,locationDir)
    cdef _propagate(particles self,double* particleArray, numpy.float64_t dt,char[] momentumDir="y", char[] locationDir="x"):
#        cdef double [:] particleArray=self.particleArray
        cdef unsigned int i = self.coordinateNames[momentumDir]
        cdef unsigned int j = self.coordinateNames[locationDir]
        cdef unsigned int totalElements=(self.particleNumber*self.dimensions)
        while j<totalElements:
            particleArray[j]+=dt*particleArray[i]
            i+=self.dimensions
            j+=self.dimensions
        
cdef class Cgrid1D:
    cdef public numpy.float64_t left
    """Minimum coordinate / left edge of the grid"""
    cdef public numpy.float64_t right
    """Maximum coordinate / right edge of the grid"""
    cdef public unsigned int nBins
    """Number of bins of the grid"""
    cdef public numpy.float64_t gridSpacing
    """Width of one grid cell"""
    cdef public numpy.ndarray density
    """Array of the densities on the grid"""
    cdef public numpy.ndarray bins
    """Array of the coordinates of the bin centers"""
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
        self.bins = numpy.arange(left,right,self.gridSpacing)+self.gridSpacing*.5
    cpdef object particles2Grid(Cgrid1D self, particles particlesIn,char[] direction="x"):
        ''' Convenience method
        allows to calculate the densities of particles on a grid.
        expects a string for the direction (direction="x", "y" or something given in the dictionary of 
        '''
        self.points2Grid(particlesIn.particleArray,particlesIn.coordinateNames[direction],pointMass=particlesIn.particleCharge)
    
    cpdef object grid2Particles(Cgrid1D self, particles particlesIn,char[] directionIn="x",char[] directionOut="y"):
        ''' Convenience method
        allows to calculate the densities of particles on a grid.
        expects a string for the direction (direction="x", "y" or something given in the dictionary of 
        '''
        self.grid2Points(particlesIn.particleArray,particlesIn.coordinateNames[directionIn],particlesIn.coordinateNames[directionOut],pointMass=particlesIn.particleCharge)
        
    cpdef object points2Grid(Cgrid1D self,numpy.ndarray[numpy.float64_t,ndim=2,mode='c'] particlesIn,direction=0,pointMass=1.):
        '''
        expects:
            a two-dimensional numpy array (in c-mode) of the form 
            [ [p1.c1, p1.c2,..., p1.cn], [p1.c2,p1.c2,..., p2.cn]... [pm.c1,pm.c2,..., pm.cn]]
            If your array is not in C-Mode (for example if you transposed it shortly before) you will likely get an error. In that case you can just
            pass a copy of your original array.
        '''
        self._points2Grid(&particlesIn[0,0],direction,particlesIn.shape[0],particlesIn.shape[1],pointMass)
    
    cpdef object grid2Points(Cgrid1D self,numpy.ndarray[numpy.float64_t,ndim=2,mode='c'] particlesIn,directionIn=0, directionOut=1,pointMass=1.):
        '''
        expects:
            a two-dimensional numpy array (in c-mode) of the form 
            [ [p1.c1, p1.c2,..., p1.cn], [p1.c2,p1.c2,..., p2.cn]... [pm.c1,pm.c2,..., pm.cn]]
            If your array is not in C-Mode (for example if you transposed it shortly before) you will likely get an error. In that case you can just
            pass a copy of your original array.
        '''
        self._grid2Points(&particlesIn[0,0], directionIn, directionOut,particlesIn.shape[0],particlesIn.shape[1],pointMass)
    
    cdef void _points2Grid(Cgrid1D self,double* particles, unsigned int direction, unsigned int numberOfParticles, unsigned int dimensions,double pointMass):
        '''
        expects:
            particles: array of coordinates (1-dimensional)
            grid: grid object
            particle coordinates should not exceed grid area
        '''
        cdef numpy.ndarray[numpy.float64_t,ndim=1] tmpgrid=numpy.zeros(self.nBins,dtype=numpy.float64)
    #    cdef double [:] particles=particlesIn
        cdef int leftIndex =0
        cdef numpy.float64_t binPosition= 0.
        cdef numpy.float64_t invBinWidth=pointMass/self.gridSpacing
        cdef numpy.float64_t left=self.left+self.gridSpacing*.5
        cdef numpy.float64_t right=self.right
        cdef numpy.float64_t gridSpacing=self.gridSpacing
        cdef unsigned int nBins=self.nBins
        cdef unsigned int i = direction

        #for i in range(2*numberOfParticles,):
        # divmod is a python function, so we are faster without it.
        cdef unsigned int totalElements=(numberOfParticles*dimensions)
        while i<totalElements:
            leftIndex=<int>((particles[i]-left)//gridSpacing) #this is the div part of divmod
            binPosition=abs((particles[i]-left)%gridSpacing)*invBinWidth#the mod part of divmod
            tmpgrid[<unsigned int>(leftIndex % nBins)]+=pointMass-binPosition #interpolation on the neighbouring bins
            tmpgrid[<unsigned int>((leftIndex + 1) % nBins)]+=binPosition
            i+=dimensions
        self.density=tmpgrid*1/(self.gridSpacing) # so far we just summed the particles, dividing by gridSpacing gives a density

    cdef void _grid2Points(Cgrid1D self,double* particles, unsigned int directionIn, unsigned int directionOut, unsigned int numberOfParticles, unsigned int dimensions,double pointMass):
        '''
        expects:
            particles: array of coordinates (1-dimensional)
            grid: grid object
            particle coordinates should not exceed grid area
        '''
        cdef numpy.ndarray[numpy.float64_t,ndim=1] tmpgrid=self.density
    #    cdef double [:] particles=particlesIn
        cdef unsigned int leftIndex =0
        cdef numpy.float64_t binPosition= 0.
        cdef numpy.float64_t invBinWidth=pointMass/self.gridSpacing
        cdef numpy.float64_t left=self.left+self.gridSpacing*.5
        cdef numpy.float64_t right=self.right
        cdef numpy.float64_t gridSpacing=self.gridSpacing
        cdef unsigned int nBins=self.nBins
        cdef unsigned int i = directionIn
        cdef unsigned int j = directionOut

        #for i in range(2*numberOfParticles,):
        # divmod is a python function, so we are faster without it.
        cdef unsigned int totalElements=(numberOfParticles*dimensions)
        while i<totalElements:
            leftIndex=<unsigned int>((particles[i]-left)//gridSpacing) #this is the div part of divmod
            binPosition=abs((particles[i]-left)%gridSpacing)*invBinWidth #the mod part of divmod
            particles[j]+=tmpgrid[(leftIndex % nBins)]*(pointMass-binPosition) #interpolation on the neighbouring bins
            particles[j]+=tmpgrid[((leftIndex + 1) % nBins)]*binPosition
            i+=dimensions
            j+=dimensions
      #      print "%d %d %d" % (i,j,(leftIndex % nBins),tmpgrid[((leftIndex + 1) % nBins)],tmpgrid[(leftIndex % nBins)])
#        self.density=tmpgrid*1/(self.gridSpacing) # so far we just summed the particles, dividing by gridSpacing gives a density
cdef class Cgrid2D:
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
    cpdef object particles2Grid(Cgrid2D self, particles pics, axisX="x",axisY="y"):
        self.points2Grid(pics.particleArray,pointMass=pics.particleCharge,direction1=pics.coordinateNames[axisX],direction2=pics.coordinateNames[axisY])
    cpdef object points2Grid(Cgrid2D self,numpy.ndarray[numpy.float64_t,ndim=2,mode='c'] particlesIn,double pointMass=1.,unsigned int direction1=0,unsigned int direction2=1):
        cdef double[:,:] density
        self.density=numpy.zeros((self.nBinsX,self.nBinsY),dtype=numpy.float64)
        density=self.density
        self._points2Grid(&particlesIn[0,0],&density[0,0],pointMass,particlesIn.shape[1],particlesIn.shape[0],direction1,direction2)
    cdef void _points2Grid(Cgrid2D self,double* particles, double* grid, double pointMass,unsigned int dimensions,unsigned int particleNumber,unsigned int direction1,unsigned int direction2):
        cdef numpy.float64_t invBinWidthX=pointMass/self.gridSpacingX
        cdef numpy.float64_t invBinWidthY=pointMass/self.gridSpacingY
        cdef unsigned int leftIndexX,leftIndexY
        cdef double binPositionX,binPositionY
        cdef double particleX,particleY
        cdef unsigned long int i=0
        cdef unsigned long int totalElements=dimensions*particleNumber-dimensions+1
        pointMassSqrt=pointMass**0.5
        cdef double subtractX=self.xmin+self.gridSpacingX/2.
        cdef double subtractY=self.ymin+self.gridSpacingY/2.
        while i <= totalElements:
            particleX=particles[<unsigned int>(i+direction1)]
            particleY=particles[<unsigned int>(i+direction2)]
            leftIndexX=<unsigned int>((particleX-subtractX)//self.gridSpacingX)
            leftIndexY=<unsigned int>((particleY-subtractY)//self.gridSpacingY)
            binPositionX=abs((particleX-subtractX)%self.gridSpacingX)*invBinWidthX
            binPositionY=abs((particleY-subtractY)%self.gridSpacingY)*invBinWidthY
            grid[self.nBinsY*(leftIndexX%self.nBinsX)+leftIndexY%self.nBinsY]+=(pointMassSqrt-binPositionX)*(pointMassSqrt-binPositionY)
            grid[self.nBinsY*((leftIndexX+1)%self.nBinsX)+leftIndexY%self.nBinsY]+=(binPositionX)*(pointMassSqrt-binPositionY)
            grid[self.nBinsY*((leftIndexX+1)%self.nBinsX)+(leftIndexY+1)%self.nBinsY]+=(binPositionX)*(binPositionY)
            grid[self.nBinsY*(leftIndexX%self.nBinsX)+(leftIndexY+1)%self.nBinsY]+=(pointMassSqrt-binPositionX)*(binPositionY)
            i+=dimensions
#        print tmpgrid
        self.density*=1/(self.gridSpacingX*self.gridSpacingY)