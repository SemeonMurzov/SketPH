import numpy as np
from sketph.data import field
class ElasticityIdealPlasticityHugoniot:

	def __init__(self, dens0, ca, sa, 
				poisonRatio, shearModulus, yieldStrength, tensileStrength,
				gamma=0, heatCapacity=1):
		self.gamma = gamma
		self.dens0 = dens0
		self.ca    = ca 
		self.sa    = sa
		self.G     = shearModulus
		self.Y     = yieldStrength
		self.tenS  = tensileStrength
		self.nu    = poisonRatio
		self.B     = dens0*ca**2
		self.Cv    = heatCapacity

	def all(self, particles):
		# bulk properties eos
		density = particles[field.density]
		energy  = particles[field.energy]
		X = self.dens0/density
		pressure = np.zeros(particles.lenActive)
		referencEnergy = np.zeros(particles.lenActive)
		ic = X<1.0
		if ic.any():
			pressure[ic] += self.B * (1.0 - X[ic]) / \
			 						(1.0 - (1.0 - X[ic])*self.sa)**2
			pressure[~ic] += self.B * (1.0 - X[~ic])
			referencEnergy[ic] = 0.5*pressure[ic]*(1.0/self.dens0 - 1.0/density[ic]) 
			referencEnergy[~ic] = 0.5*self.B*(1.0 - X[~ic])*(1.0 - X[~ic])/self.dens0
		else:
			pressure += self.B * (1.0 - X)
			referencEnergy += 0.5*self.B*(1.0 - X)*(1.0 - X)/self.dens0

		pressure+= self.gamma*density*(energy-referencEnergy)
		# semi-acoustic approximation
		cb = np.sqrt(self.ca*self.ca +
		 		4.0*self.sa*pressure/self.dens0)
		ic = cb<0.01*self.ca
		if np.any(ic):
			cb[ic] = 0.01*self.ca

		# elasticity begins 
		G =  1.5 * density * cb * cb * (1.0 - 2.0*self.nu)/(1.0 + self.nu)
		deviatorStress = particles[field.deviatorStress]
		# deviatorStressZZ = - deviatorStress[:,0] - deviatorStress[:,3]
		xxmyy = deviatorStress[:,0] - deviatorStress[:,3]
		yymzz = deviatorStress[:,0]
		zzmxx = - 2*deviatorStress[:,0] - deviatorStress[:,3]

		J2 = 0.5 * (xxmyy*xxmyy + yymzz*yymzz + zzmxx*zzmxx)\
		          + 3.0 * deviatorStress[:,1]**2

		J2 = np.sqrt(J2)
		particles[field.elpli] = False
		iplastic = J2 >= self.Y
		if iplastic.any():
			(particles[field.elpli])[iplastic] = True
			(particles[field.deviatorStress][:,0])[iplastic] *= self.Y/J2[iplastic]
			(particles[field.deviatorStress][:,1])[iplastic] *= self.Y/J2[iplastic]
			(particles[field.deviatorStress][:,2])[iplastic] *= self.Y/J2[iplastic]
			(particles[field.deviatorStress][:,3])[iplastic] *= self.Y/J2[iplastic]

		particles[field.soundSpeedL] = np.sqrt(cb*cb + 4.0*G/(3.0*density))
		particles[field.soundSpeedT] = np.sqrt(G/density)
		particles[field.soundSpeed]  = cb
		particles[field.pressure]    = pressure
		particles[field.shearModulus]= G

