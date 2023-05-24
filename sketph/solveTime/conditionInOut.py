from sketph.data import field 
import numpy as np
from sketph.box.rectangle import Box
from .euler import TimeStepper



class stepperAMW:
	def __init__(self, coreStepper=TimeStepper()):
		self.stepper = coreStepper
		pass

	def collectParticlesInOutFix(self,particles, particlesInOut, iInOut):
		NInOut = len(iInOut[0])
		particlesInOut[field.coords][:NInOut,:] = particles[field.coords][iInOut,:]
		particlesInOut[field.density][:NInOut] = particles[field.density][iInOut]
		particlesInOut[field.pressure][:NInOut] = particles[field.pressure][iInOut]
		particlesInOut[field.energy][:NInOut] = particles[field.energy][iInOut]
		particlesInOut[field.deviatorStress][:NInOut,:] = particles[field.deviatorStress][iInOut]
		particlesInOut[field.soundSpeed][:NInOut] = particles[field.soundSpeed][iInOut]
		particlesInOut[field.soundSpeedL][:NInOut] = particles[field.soundSpeedL][iInOut]
		particlesInOut[field.soundSpeedT][:NInOut] = particles[field.soundSpeedL][iInOut]
		particlesInOut[field.size][:NInOut] = particles[field.size][iInOut]

	def replaceParticlesMWInOutFix(self, particles, particlesInOut, iInOut, VInOut, dt):
		NInOut = len(iInOut[0])
		particles[field.coords][iInOut,0] = particlesInOut[field.coords][:NInOut,0]+VInOut*dt
		particles[field.coords][iInOut,1] = particlesInOut[field.coords][:NInOut,1]
		particles[field.velocity][iInOut,0] = VInOut
		particles[field.velocity][iInOut,1] = 0.0
		particles[field.density][iInOut] = particlesInOut[field.density][:NInOut]
		particles[field.pressure][iInOut] = particlesInOut[field.pressure][:NInOut]
		particles[field.energy][iInOut] = particlesInOut[field.energy][:NInOut]
		particles[field.deviatorStress][iInOut,:] = particlesInOut[field.deviatorStress][:NInOut,:]
		particles[field.soundSpeed][iInOut] = particlesInOut[field.soundSpeed][:NInOut]
		particles[field.soundSpeedL][iInOut] = particlesInOut[field.soundSpeedL][:NInOut]
		particles[field.soundSpeedT][iInOut] = particlesInOut[field.soundSpeedT][:NInOut]
		particles[field.size][iInOut] = particlesInOut[field.size][:NInOut]

	def collectParticlesInOut(self, particles, particlesInOut, iInOut):
		NInOut = len(iInOut[0])
		particlesInOut[field.coords][:NInOut,:] = particles[field.coords][iInOut,:]
		particlesInOut[field.deviatorStress][:NInOut,:] = particles[field.deviatorStress][iInOut,:]

	def replaceParticlesMWInOut(self, particles, particlesInOut, iInOut, VInOut, dt):
		NInOut = len(iInOut[0])
		particles[field.coords][iInOut,0] = particlesInOut[field.coords][:NInOut,0]+VInOut*dt
		particles[field.coords][iInOut,1] = particlesInOut[field.coords][:NInOut,1]
		particles[field.deviatorStress][iInOut,:] = particlesInOut[field.deviatorStress][:NInOut,:]
		particles[field.velocity][iInOut,0] = VInOut
		particles[field.velocity][iInOut,1] = 0.0

	def replaceParticlesMWInOutFree(self, particles, particlesInOut, iInOut, VInOut, dt):
		NInOut = len(iInOut[0])
		particles[field.coords][iInOut,0] = particlesInOut[field.coords][:NInOut,0]+VInOut*dt
		particles[field.coords][iInOut,1] = particlesInOut[field.coords][:NInOut,1]
		particles[field.velocity][iInOut,0] = VInOut
		particles[field.velocity][iInOut,1] = 0.0



	def __call__(self, particles, particlesIn, particlesOut, particlesFreeze,
								  xSliceIn,    xSliceOut,    xSliceFreeze, 
								  upIn, upOut,  dt, CFL):
		# разметка входящих и выходящих частиц
		# marks inflow/outflow particles
		iIn  = np.where(np.greater(particles[field.coords][:,0] - xSliceIn, 0))
		self.collectParticlesInOut(particles, particlesIn, iIn)
		iOut = np.where(np.less(particles[field.coords][:,0] - xSliceOut, 0))
		self.collectParticlesInOut(particles, particlesOut, iOut)
		iFreeze = np.where(np.less(particles[field.coords][:,0] - xSliceFreeze, 0))
		self.collectParticlesInOutFix(particles, particlesFreeze, iFreeze)
		# совершение шага по времени стандартным решателем
		# perform a contact SPH timestep
		self.stepper.makeStep(particles, dt)
		# корректировка положений частиц находящихся в граничной области
		# replace the boundary layered particles
		self.replaceParticlesMWInOutFree(particles, particlesIn, iIn, upIn, dt)
		self.replaceParticlesMWInOut(particles, particlesOut, iOut, upOut, dt)
		self.replaceParticlesMWInOutFix(particles, particlesFreeze, iFreeze, upOut, dt)

		