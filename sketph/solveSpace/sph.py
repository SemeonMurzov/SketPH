import numpy as np
from sketph.data import field
from sketph.math import WendlandC2
from sketph.box.rectangle import Box
from sketph.model import ElasticityIdealPlasticityHugoniot
dim = 2
# right-hand of equations 
class SpatialSum:
	def __init__(self, 
				box = Box(),
				boundaries = {0:"period", 1:"period"}, 
				kernel = WendlandC2(),
				materials = {0:ElasticityIdealPlasticityHugoniot},
				smoothingScale = 1.0):
		# smoothing kernel WendlandC2
		self.kernel = kernel
		# a list of boundary conditions
		self.boundaries = boundaries
		# simulation box parameters
		self.box = box
		""" factor of smoothing scale length 
		which detemines the number of neighbours to be included
		"""
		self.ss = smoothingScale
		self.materials = materials
		

	def __call__(self, particles, spatialNeibs):
		particles[field.volumeRate][:]      = 0
		particles[field.deviatorSRate][:,:] = 0
		particles[field.energyRate][:]      = 0
		particles[field.force][:,:]         = 0
		particles[field.angularRate][:]     = 0
		# kernel derivative
		dW = self.kernel.derivative

		e_ji = np.zeros(dim, dtype="f8")
		for i in range(particles.lenActive):
			dens_i   = particles[field.density][i]
			cl_i     = particles[field.soundSpeedL][i]
			ct_i     = particles[field.soundSpeedT][i]
			stress_i = particles[field.deviatorStress][i,:].reshape(dim,dim)-\
					   particles[field.pressure][i]*np.eye(dim)
			z_i      = np.array([dens_i*cl_i, dens_i*ct_i])
			for j in spatialNeibs.getIt(particles.data[field.coords][i]):
				if i==j:
					continue
				for k in range(dim):
					e_ji[k] = particles.data[field.coords][j,k]-\
							  particles.data[field.coords][i,k]
					if self.boundaries[k]=="period":
						if e_ji[k]>self.box.Lmax[k]:
							e_ji[k]-=self.box.L[k]
						if e_ji[k]<self.box.Lmin[k]:
							e_ji[k]+=self.box.L[k]
				r = np.sqrt(np.dot(e_ji,e_ji))

				size_i = particles[field.size][i]
				size_j = particles[field.size][j]
				h = self.ss*(size_i+size_j)
				if r>h:
					continue
				e_ji = e_ji/r
				rotator =  np.array([[ e_ji[0], e_ji[1]],
									 [-e_ji[1], e_ji[0]]], dtype="f8")
				u_i      = rotator@particles[field.velocity][i]
				u_j      = rotator@particles[field.velocity][j]
				dens_j   = particles[field.density][j]
				mass_j   = particles[field.mass][j]
				cl_j     = particles[field.soundSpeedL][j]
				ct_j     = particles[field.soundSpeedT][j]
				stress_j = particles[field.deviatorStress][j,:].reshape(dim,dim)-\
						   particles[field.pressure][j]*np.eye(dim)
				z_j      = np.array([dens_j*cl_j, dens_j*ct_j])
				
				projStress_i = (rotator@stress_i)@e_ji
				projStress_j = (rotator@stress_j)@e_ji
				# to particle i from particle j
				dWji = dW(r, h) * mass_j / dens_j

				# vector of stress tensor  components projected to the the radial line and it's |_ line sums 
				stressContact = (projStress_j * z_i + projStress_i * z_j + 
								(u_j - u_i) * z_i * z_j) / (z_i + z_j)	
				# correspondent velocity components
				uContact = (u_i * z_i + u_j * z_j + (projStress_j - projStress_i)) / (z_i + z_j)



				# check adhesion/friction at materials interface

				reverseRotator = np.transpose(rotator)
				u_Ri = uContact - u_i
				v_Ri = reverseRotator @ u_Ri
				particles[field.force][i]      -= dWji * (reverseRotator @ stressContact)
				particles[field.volumeRate][i] -= dWji * u_Ri[0]
				particles[field.energyRate][i] -= dWji * np.dot(stressContact, uContact)
				# elasticity deviatoric strain rate accumulation
				particles[field.deviatorSRate][i,0] -= dWji * 2.0 * (3.0 * v_Ri[0] * e_ji[0] - u_Ri[0]) / 3.0
				particles[field.deviatorSRate][i,1] -= dWji * (v_Ri[0] * e_ji[1] + v_Ri[1] * e_ji[0])
				particles[field.deviatorSRate][i,3] -= dWji * 2.0 * (3.0 * v_Ri[1] * e_ji[1] - u_Ri[0]) / 3.0

				# local angular velocity estimate 
				# to make a correction to the rigid rotation of the deviatoric stress tensor
				particles[field.angularRate][i]   -= dWji*(v_Ri[1] * e_ji[0] - v_Ri[0] * e_ji[1])
		# convert to acceleration
		particles[field.force][:,0] *= 2.0 / particles[field.density][:]
		particles[field.force][:,1] *= 2.0 / particles[field.density][:]

		particles[field.volumeRate][:] *= 2.0
		# convert to energy density rate
		particles[field.energyRate][:] *= 2.0 / particles[field.density][:]

		# convert to deviatoric stress rate
		particles[field.deviatorSRate][:,0] *= 2.0*particles[field.shearModulus][:]
		particles[field.deviatorSRate][:,1] *= 2.0*particles[field.shearModulus][:]
		particles[field.deviatorSRate][:,2]  = particles[field.deviatorSRate][:,1]
		particles[field.deviatorSRate][:,3] *= 2.0*particles[field.shearModulus][:]
		#TODO sums to j in irder to make a loop over i>j only
		