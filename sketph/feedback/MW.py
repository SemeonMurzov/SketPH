import numpy as np
import collections 
import statistics as stat
import os
from sketph.box.rectangle import Box
from sketph.data import field

dim = 2
def leastSquareL2Residue(x,y,LSCoefs):
	return np.sum((y-LSCoefs[0]*x-LSCoefs[1])**2)
    
def leastSquareResidue(x,y,LSCoefs):
	return y-LSCoefs[0]*x-LSCoefs[1]

class MW:
	def __init__(self, box = Box(), dirname = "outdat", nStartBuffer = 25, 
				  frontPositionCoef = 0.8, target = 1, sigmaWish = 3.0e8, up=1,
				  kappa = 4250, delta = 10):
		# среднее напряжение и значение целевой функции в поиске
		# volume average of XX stress norm as a target value of the target
		self.sigmaWish = sigmaWish #Pa
		self.target = target
		self.delta = delta
		self.kappa = kappa
		self.up = up
        # длина памяти буфера начальная
		# initial memory length of the buffer
		self.nBuffer = nStartBuffer
		# длина буфера
		# буферы с историей
		# history buffers
		self.Buffers = {}
		self.BuffersNames = ["targetsB", "timesB", "upOut"]
		for Name in self.BuffersNames:
			self.Buffers[Name] = collections.deque(maxlen = self.nBuffer)
		if not os.path.exists(dirname):
			os.mkdir(dirname)
		self.logfile = open(dirname+"//log.dat","w")


	# less iffness and selfness
	def __call__(self, particles, time, 
					  upIn, upOut, dt, nSteps):
		mechanicalEnergyAverage = np.sum((particles[field.pressure]-
										particles[field.deviatorStress][:,0])*
										particles[field.size]**dim)
		volume  = np.sum(particles[field.size]**dim)
		sigmaAverage = mechanicalEnergyAverage/volume
		omega = np.max([0.0,sigmaAverage/self.sigmaWish])
		self.Buffers["timesB"].append(time)
		self.Buffers["targetsB"].append(omega)
		bell  = 1.0/(1.0+(self.delta*abs(omega-self.target)))

		omega = np.mean(np.array(self.Buffers["targetsB"])[-self.nBuffer:])
		upOut = max(upIn, min(-self.kappa*(1.0+omega-(omega-self.target)* \
		                             bell)/2.0,-self.kappa/1.3) )
		messageOut = str(omega)+"\t"+str(time)+"\t"+str(upOut)+"\t"+str(nSteps)
		self.logfile.write(messageOut+"\n")
		print(messageOut)

		self.Buffers["upOut"].append(upOut)
		return upOut

	
	def closeLog(self):
		self.logfile.close()