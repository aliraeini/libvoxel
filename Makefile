
srcs := voxelImage.cpp vxlPro.cpp
mains := Ufraw2Uc  voxelToFoamPar  voxelToFoam  voxelImageProcess
all: $(mains)
tsts := test.py
$(info mains: $(mains),   tsts: $(tsts) )


#ifneq ("${OPT}", ".exe")
#  msCFLAGS += -DLPNG
#  msLFLAGS += -lpng
#endif
msSrc ?= $(abspath ..)
USE_ZLIB=1
USE_TIFF=1
#USE_SVG=1
USE_OMP=1
MS_VERSION=2
USE_SINGLECPP=1
USE_msTEST=1
include  ${msSrc}/script/Makefile.in
