import sys
import os
# data
# classes
from sketph.box.rectangle import Box
from sketph.box.rectangle import addPlate
from sketph.model import ElasticityIdealPlasticityHugoniot


def master_print(s):
	print(s)
	sys.stdout.flush()
def printComplete(s):
	print(s+"...")

dirname = "inlet"
if (len(sys.argv) > 1):
	dirname = sys.argv[1]
if not os.path.exists(dirname):
	os.mkdir(dirname)

# material parameters
copper_dens0 = 8960.
copper_ca = 3930.
bulk_modulus = copper_dens0*copper_ca**2
shear_modulus = 4.3e10
copper_poissonRatio = (3*bulk_modulus-2*shear_modulus)/(2*(3*bulk_modulus+shear_modulus))
a_copper = 1.5
gamma_copper = 1.7
heat_capacity = 384.603
yield_strength_copper = 3.5e+8
tensile_stress= 3.5e+8

material = ElasticityIdealPlasticityHugoniot(
	copper_dens0, copper_ca, a_copper, 
	copper_poissonRatio, shear_modulus,
	yield_strength_copper, tensile_stress,
	gamma=gamma_copper, heatCapacity=heat_capacity)

# simulation box dimensions
dim = 2
maxDim = 6e-6
D = 2e-6
# replicas number of initial sample 
NxBox = 1
# rectangle simulation box  
xMax = NxBox*maxDim
xMin =-NxBox*maxDim
yMin =-maxDim
yMax = maxDim
# class simulation box
box = Box(dim = 2,shape = "rectangular",
		box = {"xmin":1.4*xMin, "xmax":1.4*xMax,
			   "ymin":yMin, "ymax":yMax})

sample = Box(dim = 2,shape = "rectangular",
		box = {"xmin":xMin, "xmax":xMax,
			   "ymin":yMin, "ymax":yMax})
# particles initialization
particles = addPlate(D, sample, dens0 = copper_dens0)
material.all(particles)
import pickle
with open(dirname+"//inlet.pkl", "wb") as f:
	pickle.dump(particles.as_dict(),f)

