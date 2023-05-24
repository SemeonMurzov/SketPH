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

class AMW:
	def __init__(self, box = Box(), dirname = "outdat", nStartBuffer = 1500, 
				  frontPositionCoef = 0.8, target = 1, sigmaWish = 3.0e8, up=1):
		# среднее напряжение и значение целевой функции в поиске
		# volume average of XX stress norm as a target value of the target
		self.sigmaWish = sigmaWish #Pa
		self.target = target
		self.up = up
        # длина памяти буфера начальная
		# initial memory length of the buffer
		self.nStartBuffer = nStartBuffer
		# количество шагов до перехода в режим стационарного течения 
		# number of steps to start a stationary regime
		self.nStationaryGlobal = 2*nStartBuffer
		# длина буфера
		# actual buffer length
		self.nBuffer = self.nStationaryGlobal
		# управляющий параметр подстройки адаптивного подвижного окна
		# the minus length form the removing plane to the shock front as
		# velocity change value control parameter of the adaptive moving window
		self.kL = -frontPositionCoef*box.L[0]
		self.OmegaA = collections.deque()#maxlen=2
		self.dOmegadT = collections.deque()#maxlen=2
		self.BufferUpIn = collections.deque()#maxlen=2
		self.nStepdOmegadT = collections.deque(maxlen=2)
		# режимы адаптивного окна
		# AMW regimes 
		# активная подстройка СО
		# active search of the coordinate system velocity
		self.mimicry = False
		# подстройка положения фронта УВ
		# shock wave front position tuning
		self.positionTuning = 0
		# изменение скорости СО при подстройке положения
		# the velocity change in the position tuning control  
		self.deltaUpPositionTuning = 0
		# буферы с историей
		# history buffers
		self.Buffers = {}
		self.BuffersIn = {}
		self.BuffersOut = {}
		self.BuffersNames = ["targetsB", "timesB", "upOut"]
		self.BuffersNamesIn = ["j","Upx","dens"]
		self.BuffersNamesOut = ["j", "Upx", "dens"]
		for Name in self.BuffersNames:
			self.Buffers[Name] = collections.deque(maxlen = self.nBuffer)
		for Name in self.BuffersNamesOut:
			self.BuffersOut[Name] = collections.deque(maxlen = self.nBuffer)
		for Name in self.BuffersNamesIn:
			self.BuffersIn[Name] = collections.deque(maxlen = self.nBuffer)
		if not os.path.exists(dirname):
			os.mkdir(dirname)
		self.logfile = open(dirname+"//log.dat","w")

	def resetBuffers(self, time, omega, upOut):
		self.Buffers["timesB"].clear()
		self.Buffers["targetsB"].clear()
		self.Buffers["timesB"].append(time)
		self.Buffers["targetsB"].append(omega)
		self.Buffers["upOut"].clear()
		self.Buffers["upOut"].append(upOut)

	# iffness and selfness
	def __call__(self, particles, time, 
					  upIn, upOut, dt, nSteps):
		estimated_model = [0,0]
		deltaUpGlobal = 0
		mechanicalEnergyAverage = np.sum((particles[field.pressure]-
										particles[field.deviatorStress][:,0])*
										particles[field.size]**dim)
		volume  = np.sum(particles[field.size]**dim)
		sigmaAverage = mechanicalEnergyAverage/volume
		omega = np.max([0.0,sigmaAverage/self.sigmaWish])
		mechanicalEnergyAverage = np.sum(particles[field.pressure]*particles[field.size]**dim)
		presAverage = mechanicalEnergyAverage/volume
		omega1 = np.max([0.0,presAverage/self.sigmaWish])
		omega2 = np.max([0.0,particles.lenActive])
		if omega>0.8*self.target or self.mimicry:
			self.Buffers["timesB"].append(time)
			self.Buffers["targetsB"].append(omega)
			print(len(self.Buffers["targetsB"]), " of ", self.nStartBuffer)
		if self.mimicry:
			# make a linear regression fit of averaged measured target function
			estimated_model = np.polyfit(self.Buffers["timesB"],self.Buffers["targetsB"],1)
			# current number of steps after the last estimate of model parameters
			self.nStationary = nSteps-self.nStepdOmegadT[-1]
			currentResidue = 0.5*np.max(np.abs(leastSquareResidue(
												np.array(self.Buffers["timesB"]),
                                                np.array(self.Buffers["targetsB"]),
												estimated_model)  )  ) 
			currentDrift = abs(estimated_model[0]*(self.Buffers["timesB"][-1]-
												   self.Buffers["timesB"][0]))
			if (currentResidue<currentDrift and self.nStationary>self.nStartBuffer) :
				self.OmegaA.append(abs(self.target-abs(self.target-stat.mean(
												   self.Buffers["targetsB"]))))
				self.dOmegadT.append(estimated_model[0])
				self.BufferUpIn.append(upIn)
				self.nStepdOmegadT.append(nSteps)
				self.resetBuffers(time, omega, upOut)
				if abs(omega-self.target)>1e-2 and abs(self.dOmegadT[-1]*self.kL)>100:
					deltaUpAdditional = (omega-self.target)*self.up
				else:
					deltaUpAdditional = 0
				deltaUpGlobal = self.dOmegadT[-1]*self.kL+deltaUpAdditional
			else:
				if self.nStationary>self.nStartBuffer:
					print(self.nStationaryGlobal-self.nStationary, 
						  " steps to stationary mode with nBuffer = ", self.nBuffer)
					print("drift<residue ",  currentDrift," < ",currentResidue)
				else:
					print("drift ", currentDrift, "\t residue ", currentResidue)
				# increase buffer length 
				if self.nStationary>=self.nBuffer:
					self.nBuffer += self.nStartBuffer
					self.Buffers["timesB"] = collections.deque(self.Buffers["timesB"],
											      maxlen = self.nBuffer)
					self.Buffers["targetsB"] = collections.deque(self.Buffers["targetsB"],
												  maxlen = self.nBuffer)
					print("buffer increased", self.nBuffer)
				# turn on the stationarity regime
				if self.nStationary == self.nStationaryGlobal:
					print("stationarity by velocity: "+str(upIn))
					self.nStationaryGlobal += self.nStartBuffer
					self.nStartBuffer = int(1.2*self.nStartBuffer)
					# additional tuning of position for even correct velocity
					self.deltaUpPositionTuning = np.sign(self.target-omega)*\
										max(abs(10*(self.target-omega)*self.up),100)
					if self.deltaUpPositionTuning+upOut>0:
						self.deltaUpPositionTuning=-upOut
					if abs(self.Buffers["targetsB"][-1]-self.target)>1e-2:
						deltaUpGlobal += self.deltaUpPositionTuning
						self.positionTuning = True
						self.mimicry = False
						self.Buffers = self.resetBuffers(time, omega, upOut)
						self.nStepdOmegadT.append(nSteps)
						print("into position tuning")
		else:
			if self.positionTuning and abs(omega-self.target)<1e-3:
				deltaUpGlobal -= self.deltaUpPositionTuning
				self.dOmegadT.append(np.polyfit(self.Buffers["timesB"],
											self.Buffers["targetsB"],1)[0])
				self.BufferUpIn.append(upIn)
				self.nStepdOmegadT.append(nSteps)
				self.resetBuffers(time, omega, upOut)
				self.positionTuning = False
				self.mimicry = True
				print("back from position tuning to mimicring of velocity")
			else:
            	#initialization after start point of algorithm 
				if omega>=self.target and (not self.positionTuning):
					self.dOmegadT.append(np.polyfit(self.Buffers["timesB"],
											 self.Buffers["targetsB"],1)[0])
					self.BufferUpIn.append(upIn)
					self.OmegaA.append(abs(1-abs(1-stat.mean(self.Buffers["targetsB"]))))
					deltaUpGlobal = self.dOmegadT[-1]*self.kL
					self.nStepdOmegadT.append(nSteps)
					self.resetBuffers(time, omega, upOut)
					self.mimicry = True
					print(self.kL*self.dOmegadT[-1], ' mimicry is started')
		if self.dOmegadT:
			a = str(self.dOmegadT[-1])
		else:
			a = "0"
		if self.mimicry:
			typeOfTuning = "m"
		else:
			if self.positionTuning!=0:
				typeOfTuning = "p"+str(self.positionTuning)
			else:
				typeOfTuning = "s"
		if estimated_model[0]:
			currentRate = estimated_model[0]
		else:
			currentRate = 0 
		if estimated_model[1]:
			currentShift = estimated_model[1]
		else:
			currentShift = 0 
		self.logfile.write(str(omega)+"\t"+str(omega1)+"\t"+str(omega2)+"\t"+
					   str(time)+"\t"+str(upIn)+"\t"+
					   str(nSteps)+"\t"+str(a)+"\t"+str(deltaUpGlobal)+"\t"+
					   str(currentRate)+"\t"+str(currentShift)+"\t"+
					   str(self.nBuffer)+"\t"+str(typeOfTuning)+"\n")
		if nSteps%10==0:
			print(str(omega)+"\t"+str(time)+"\t"+str(upIn)+"\t"+
			str(nSteps)+"\t"+str(a)+"\t"+str(deltaUpGlobal)+"\t"+str(typeOfTuning))
		if deltaUpGlobal>-upOut:
			print("unstable estimate: ", upOut, deltaUpGlobal, deltaUpAdditional)
			deltaUpGlobal = min(min(deltaUpGlobal,self.dOmegadT[-1]),-upOut-300)
			print("deltaUp ", deltaUpGlobal)
		self.Buffers["upOut"].append(upOut)
		return deltaUpGlobal
	
	def closeLog(self):
		self.logfile.close()