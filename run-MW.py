import numpy as np
import sys
import os
import matplotlib.pyplot as plt 
import pickle

# data
from sketph.data import field
from sketph.data import Storage
# algorithms
from sketph.math import WendlandC2
from sketph.neibs.spatialNeibs import SpatialNeibs
from sketph.box.rectangle import Box
from sketph.box.rectangle import replicatePlate
from sketph.model import ElasticityIdealPlasticityHugoniot
from sketph.solveSpace.sph import SpatialSum
from sketph.solveTime.euler import TimeStepper
from sketph.solveTime.conditionInOut import stepperAMW
from sketph.io.png.plotter import Plot
from sketph.feed import insertSorted
from sketph.feedback.MW import MW
"""
A sequence of classes the computational algorithm:
contains simulation box parameters and it's boundary conditions type.
the neighbors list builds according to the simulation box dimensions and
the SPH-sum of spatial derivatives to obtain time-derivatives for time step 
All of these are merged into the one solver in order to avoid
if-else "unified" approach, but it's done for a good reason - to make the code 
more explicit as a full problem, but also as a constructor of modules.
Even if the "bad" global variables are contained in sketph.data module, where 
SPH-particles field. If you want to make it global make it small name hidden in the 
class initialization definition(no, but yes).
The Storage class is defined in the same module to manage data.
To sum up, a readable minimal program in order to check computational method. 
"""

dirname = "results-png"
if (len(sys.argv) > 1):
	dirname = sys.argv[1]
if not os.path.exists(dirname):
	os.mkdir(dirname)
dirname1 = "results-pkl"
if not os.path.exists(dirname1):
	os.mkdir(dirname1)


# material parameters
copper_dens0 = 8960.
copper_ca = 3930.
bulk_modulus = copper_dens0*copper_ca**2
shear_modulus = 4.3e10
copper_poissonRatio = (3*bulk_modulus-2*shear_modulus)/\
					  (2*(3*bulk_modulus+shear_modulus))
a_copper = 1.5
gamma_copper = 1.7
heat_capacity = 384.603
yield_strength_copper = 3.5e+8
tensile_stress= 3.5e+8

materialClosure = ElasticityIdealPlasticityHugoniot(
	copper_dens0, copper_ca, a_copper, 
	copper_poissonRatio, shear_modulus,
	yield_strength_copper, tensile_stress,
	gamma=gamma_copper, heatCapacity=heat_capacity)

# simulation box dimensions
dim = 2
maxDim = 6e-6
D = 2e-6
# replica's of initial sample number 
NxBox = 18
# rectangle simulation box  
xMax = NxBox*maxDim
xMin =-NxBox*maxDim
yMin =-maxDim
yMax = maxDim

# minus velocity after the shock in laboratory system
# velocity change in a shock wave
up = -1000
# outflow velocity
upOut = -4250
# inflow velocity
upIn = up+upOut
# global (all) particles velocity change
deltaUpGlobal = 0
sigmaWish = 25e9
frontPositionCoef  = 0.5
nStartBuffer = 200
# class simulation box
box = Box(dim = 2,shape = "rectangular",
		box = {"xmin":1.4*xMin, "xmax":1.4*xMax,
			   "ymin":yMin, "ymax":yMax})

sample = Box(dim = 2,shape = "rectangular",
		box = {"xmin":xMin, "xmax":xMax,
			   "ymin":yMin, "ymax":yMax})
# particles initialization
import pickle
with open("inlet//inlet.pkl",'rb') as inFile:
	inflowSample = Storage()
	inflowSample.load(pickle.load(inFile))

particles = replicatePlate(sample, inflowSample, NxBox, maxDim, copper_dens0)
# tranform coordssetting left boundary X=0 in inflow sample for inserting
inflowSample[field.coords][:,0] += maxDim
particlesIn  = Storage(field.fields_list_in_out, int(0.3*particles.len))
print("Inital particles number ", particlesIn.lenActive)
particlesOut = Storage(field.fields_list_in_out, int(0.3*particles.len))
particlesFreeze = Storage(field.fields_list_in_out, int(0.2*particles.len))
print("particles allocated")
 
materialClosure.all(particles)

# initial conditions
# TODO restart if start_new:
particles[field.velocity][:,0] = upIn
particles[field.velocity][:,1] = 0

# inflow boundary
inflowPosition = xMax

# if start_new:
xInsert = xMax
sliceIndex = 0
xLast = 0


# inflow plane position
xMWOut = xMin
xMWIn  = xMax
# character buffer zone width
kD = 2*D
# positions of inflow/outflow zone boundaries 
xSliceIn = xMax - 2*kD
xSliceOut= xMin + 2*kD
xSliceFreeze = xSliceOut-1*kD

# dictionary passes to classes where these boundaries used
boundaries = {0:"period", 1:"period"}
neibsTable = SpatialNeibs(box = box,
			particles  = particles, 
			smoothingScale = 1.,
			buf = 2.5)

# contact smoothed particle hydrodynamics
# constructor spatial solver
spaceSum = SpatialSum(
			box = box,
			boundaries = boundaries, 
			kernel = WendlandC2(),
			materials = {0:materialClosure},
			smoothingScale = 1.0)

# euler time step
SPHstepper = TimeStepper(
			box=box, 
			boundaries = boundaries)

# pre- and post- step corrections
stepper = stepperAMW(coreStepper = SPHstepper)

# moving window constructor
MovingWindow = MW(box = sample, dirname = "logAMW", 
								nStartBuffer = nStartBuffer,
				   frontPositionCoef = frontPositionCoef, up=up,
				   			  target = 1, sigmaWish = sigmaWish)
time = 0.0
nSteps = 0
nstepsOut = 10
endStep = 2000+nSteps
nstepUpdate = 2

scatter = False
plots = True
if scatter:
	plotScatter = Plot(particles, figsize=(np.min([NxBox*4,75]),5))
	plotScatter.particles(
		xmin = xMin,
		xmax = xMax,
		ymin = yMin,
		ymax = yMax,
		timestep = nSteps
	)
if not scatter and plots:
	plot = Plot(particles, figsize=(12,8))

while nSteps < endStep:

	if nSteps<100:
		CFL = 0.25 
	else:
		CFL = 0.5
	dt = SPHstepper.getTimeStep(particles)

	spaceSum(particles,neibsTable)
	stepper(particles, particlesIn, particlesOut, particlesFreeze, 
					   xSliceIn,    xSliceOut,    xSliceFreeze, 
					   upIn, upOut,  dt, CFL)
	xInsert += upIn*dt
	
	if nSteps%nstepUpdate == 0:
		if nSteps>0:
			xLast, sliceIndex, xInsert = insertSorted.feed(particles, neibsTable, upIn, xInsert, 
						 inflowSample=inflowSample, firstInsertIndex=sliceIndex, xStart = xLast,
						 boxInflowPeriod=2*maxDim, inflowWindowPosition=inflowPosition)
	upOut = MovingWindow(particles, time, upIn, upOut, dt, nSteps)
	idelete = np.where(particles[field.coords][:,0]<xMWOut)
	if idelete[0].size>0:
		# return the indexes of particle which are passive now
		ideleted = particles.delete(idelete[0])
		''' remove the links of that nodes to the other elements of nodes
			in order to isolate these particles forcely, otherwisely
			although these particels are not still in particles slice return
			there are still nodes with indexes of particles, which removed from
			from the calculation and stays in an additional space after Storage.lenActive
		'''
		neibsTable.replaceNodes(ideleted,idelete[0])
	if nSteps%nstepUpdate == 0:
		neibsTable.checkUpdate(particles)


	time  += dt
	nSteps+= 1
	if nSteps%(5*nstepsOut)==0:
		with open("results-pkl//save_"+str(nSteps)+".pkl", "wb") as f:
			pickle.dump(particles.as_dict(),f)
		with open("results-pkl//time_"+str(nSteps)+".dat", "w") as f:
			f.write(str(time))

	if scatter:
		if nSteps%nstepsOut==0:
			plotScatter.particles(
				xmin = xMin,
				xmax = xMax,
				ymin = yMin,
				ymax = yMax,
				timestep = nSteps
			)

	if not scatter and plots:
		if nSteps%nstepsOut==0:
			for yfield, component in [(field.density,None),(field.energy, None), (field.elpli,None),
									(field.pressure,None),(field.velocity,0),(field.velocity,1),
									(field.force,0),(field.force,1),
									(field.deviatorSRate,0),(field.deviatorSRate,1),
									(field.deviatorSRate,3),(field.deviatorStress,0),
									(field.deviatorStress,1),(field.deviatorStress,3),
									(field.angularRate,None)]: 
				plot.render(
					xfield=field.coords,
					yfield=yfield,
					component = component,
					xmin = xMin,
					xmax = xMax
				)
				if component:
					dirname1 = dirname+"//"+yfield+"_"+str(component)
					if not os.path.exists(dirname1):
						os.mkdir(dirname1)
					plot.save(dirname1+"//"+yfield+"_"+str(component)+"_"+str(nSteps)+".png", fmt='png')
				else:
					dirname1 = dirname+"//"+yfield
					if not os.path.exists(dirname1):
						os.mkdir(dirname1)
					plot.save(dirname1+"//"+yfield+"_"+str(nSteps)+".png", fmt='png')
				plot.clear()
			plot.renderS(
				xfield=field.coords,
				component = 0,
				xmin = xMin,
				xmax = xMax
			)
			plot.save(dirname+"//"+"stressxx"+"_"+str(nSteps)+".png", fmt='png')
			plot.clear()


	# Hugoniot EOS,ideal ideal plasticity
	materialClosure.all(particles)
adaptiveMovingWindow.closeLog()
