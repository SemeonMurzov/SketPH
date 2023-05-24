from sketph.data import field 
import numpy as np
from sketph.box.rectangle import Box
dim = 2
class TimeStepper:
	def __init__(self, 
				box = Box(),
				boundaries = {0:"period", 1:"period"}
				):
		self.box = box
		self.boundaries = boundaries
	def getTimeStep(self, particles, 
				CFL = 0.5,):
		c2 = particles[field.soundSpeedL]**2
		s2 = (particles[field.size]*particles[field.volumeRate])**2
		il = np.where(particles[field.volumeRate]<0)
		s2[il] *= 256
		new_dt = np.amin(CFL*particles[field.size]/np.sqrt(c2+s2))
		return np.min([1e15,new_dt])

	def makeStep(self, particles, dt):
		vnew = particles[field.velocity] + particles[field.force]*dt
		particles[field.density] *= np.exp(-dt*particles[field.volumeRate])
		particles[field.energy]  += dt*particles[field.energyRate] - \
									0.5*(vnew[:,0]*vnew[:,0]+
										 vnew[:,1]*vnew[:,1]- 
										 particles[field.velocity][:,0]*
										 particles[field.velocity][:,0]-
										 particles[field.velocity][:,1]*
										 particles[field.velocity][:,1])
		particles[field.coords]  += (particles[field.velocity] + vnew)*(0.5*dt)
		
		for k in range(dim):
			if self.boundaries[k]=="period":
				im = particles[field.coords][:,k]<self.box.Lmin[k]
				particles[field.coords][im,k] += self.box.L[k]
				ip = particles[field.coords][:,k]>self.box.Lmax[k]
				particles[field.coords][ip,k] -= self.box.L[k]
		particles[field.velocity] = vnew
		stressDeviator = particles[field.deviatorStress][:]+dt*particles[field.deviatorSRate][:]
		# correct stressDeviator according to particle rotation as a rigid body

		omegadt = particles[field.angularRate]*dt
		omegadt2 = omegadt**2

		particles[field.deviatorStress][:,0]  = stressDeviator[:,0]+2*omegadt*stressDeviator[:,1]\
 																	+omegadt2*stressDeviator[:,3]
		particles[field.deviatorStress][:,1]  = stressDeviator[:,1]*(1-omegadt2)+omegadt*(
 												 stressDeviator[:,3]-stressDeviator[:,0])
		particles[field.deviatorStress][:,2]  = particles[field.deviatorStress][:,1]
		particles[field.deviatorStress][:,3]  = stressDeviator[:,3]-2*omegadt*stressDeviator[:,1]\
 																	+omegadt2*stressDeviator[:,0]
