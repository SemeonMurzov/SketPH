# ядро WendlandC2
import numpy as np


class WendlandC2:
#2D
	def __call__(self, r, h):
		hc = h*9.48683298050513768018277e-01
		q = r / hc
		if q<1:
			return 7.0/np.pi*(1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 - q)*(1.0 + 4.0*q)/hc/hc
		else:
			return 0
	def derivative(self, r, h):
		hc = h*9.48683298050513768018277e-01
		q = r / hc
		if q<1:
			return -140.0/np.pi*q*(1.0 - q)*(1.0 - q)*(1.0 - q)/hc/hc/hc
		else:
			return 0