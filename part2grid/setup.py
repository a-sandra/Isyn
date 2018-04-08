from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
# distutils: libraries = cerf
# distutils: include_dirs = /home/paulg/local/include

ext_modules = [ Extension("particles2grid",["particles2grid.pyx"],include_dirs=[numpy.get_include()],libraries=["m"], extra_compile_args=["-O4","-ftree-vectorize", "-msse2"])]#,Extension("RFcavity",["RFcavity.pyx"])]#, libraries=["cerf"], library_dirs=["/home/paulg/local/lib"],include_dirs = ['.','/home/paulg/local/include/'])]

setup(
      name='cython test class',
      cmdclass={'build_ext':build_ext},
      ext_modules=ext_modules
      )