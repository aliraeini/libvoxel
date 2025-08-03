
import os; ''' ========== set up paths  =========== '''
if not ("msRoot" in os.environ):
  print("try again after running:\nsource .../src/bashrc"); exit(-1);
from msrc  import *


nErrs=0

with open("voxcylinder.mhd", 'w') as f1:
	f1.write("""DimSize = 20 20 20
		Offset =      0    0    0
		replaceRange 0 255 1
		reset  dx 1
		Paint cylinder 0 10 10   20 10 10  5
		reset  dx 1e-6
		""");#ElementDataFile = NO_READ

runSh('.', "voxelImageProcess voxcylinder.mhd voxcylinder.tif");
runSh('.', "voxelImageProcess voxcylinder.tif voxcylinder.mhd");
totalPorosity = math.pi*5*5/(20*20)
if fileFloatDiffersFrom("voxelImageProcess.log","total_porosity:", totalPorosity, 0.05):
 if fileFloatDiffersFrom("voxelImageProcess.log","totalPorosity:", totalPorosity, 0.05): # TODO: remove cross-compatibility
    nErrs+=1

exit(nErrs)
