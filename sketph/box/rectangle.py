import numpy as np
from sketph.data import field
from sketph.data import Storage
from sketph.data.field import fields_list
import sys
class Box:
	def __init__(self,dim = 2,shape = "rectangular",
			box   = {"xmin":0,"xmax":1, "ymin":0, "ymax":1} ):
		if shape=="rectangular":
			self.Lmin = np.zeros(dim,dtype="<f8")
			self.Lmax = np.zeros(dim,dtype="<f8")
			self.L    = np.zeros(dim,dtype="<f8")
			self.Lmin[0] = box["xmin"]
			self.Lmax[0] = box["xmax"]
			if dim>1:
				self.Lmin[1] = box["ymin"]
				self.Lmax[1] = box["ymax"]
			for k in range(dim):
				self.L[k] = self.Lmax[k]-self.Lmin[k]
		else:
			raise Exception('box shape ', shape, ' should be implemented yet')


def addPlate(D, box, dens0=1):
	nX = int((box.Lmax[0] - box.Lmin[0])/D)
	nY = int((box.Lmax[1] - box.Lmin[1])/D)
	[X, Y] = np.meshgrid(
		np.linspace(box.Lmin[0] + 0.25*D, box.Lmax[0] - 0.75*D, nX),
		np.linspace(box.Lmin[1] + 0.25*D, box.Lmax[1] - 0.75*D, nY),
	)
	Nparticles =  2*len(X.ravel())
	print("number of aprticles in a sample: ", Nparticles)
	if 2*len(X.ravel())>1e7:
		print('too large')
		sys.exit()
	particlesToAdd = Storage(fields_list, Nparticles)
	unsorted = np.zeros(2*len(X.ravel()))
	unsorted[::2] = X.ravel()
	unsorted[1::2] = X.ravel()+0.5*D
	sorted_indexes = unsorted.argsort()
	particlesToAdd[field.coords][:,0] = unsorted[sorted_indexes]
	unsorted[::2] = Y.ravel()
	unsorted[1::2] = Y.ravel() + 0.5*D
	particlesToAdd[field.coords][:,1] = unsorted[sorted_indexes]
	particlesToAdd[field.size] = D/np.sqrt(2)
	particlesToAdd[field.material] = 0
	particlesToAdd[field.density] = dens0
	particlesToAdd[field.mass] = dens0*particlesToAdd[field.size]**2
	particlesToAdd[field.energy] = 0
	return particlesToAdd

def replicatePlate(box, baseSample, NxBox, maxDim, dens0):
	if NxBox*baseSample.lenActive>1e7:
		print('too large')
		sys.exit()
	particles = Storage(fields_list,NxBox*baseSample.lenActive)
	for i in range(NxBox):
		particles[field.coords][i*baseSample.lenActive:baseSample.lenActive*(i+1),0] = \
															baseSample[field.coords][:,0]+\
															 2*i*maxDim+(box.Lmin[0]+maxDim)
		particles[field.coords][i*baseSample.lenActive:baseSample.lenActive*(i+1),1] = \
															baseSample[field.coords][:,1]
		particles[field.size][i*baseSample.lenActive:baseSample.lenActive*(i+1)] = \
															baseSample[field.size]
	particles[field.material] = 0
	particles[field.density] = dens0
	particles[field.mass] = dens0*particles[field.size]**2
	particles[field.energy] = 0
	return particles