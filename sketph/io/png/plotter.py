import matplotlib.pyplot as plt
from sketph.data import field

class Plot:
	
	def __init__(self, storage, title="", figsize = (12, 8)):
		self.storage = storage
		self.title = title
		self.fig, self.axes = plt.subplots(figsize=figsize)

	def save(self, filename, fmt='png'):
		plt.savefig(filename, format=fmt)
	
	def show(self):
		plt.show()

	def get_min(self, field, component=None):
		xmin = 1e+15
		if component is None:
			for i in range(self.storage.lenActive):
				xmin = min(xmin, self.storage[field][i])
		else:
			for i in range(self.storage.lenActive):
				xmin = min(xmin, self.storage[field][i,component])
	def get_max(self, field, component = None):
		xmax = -1e+15
		if component is None:
			for i in range(self.storage.lenActive):
				xmax = max(xmax, self.storage[field][i])
		else:
			for i in range(self.storage.lenActive):
				xmax = max(xmax, self.storage[field][i,component])


	def render(self, yfield = field.velocity, xfield=field.coords,
				component = None, xmin = None, xmax = None, 
							   ymin = None, ymax = None):
		if len(self.storage[yfield].shape)>1:
			self.axes.plot(self.storage[xfield][:,0], self.storage[yfield][:,component], 
				label=yfield, linestyle="", marker='.')
		else:
			self.axes.plot(self.storage[xfield][:,0], self.storage[yfield], 
				label=yfield, linestyle="", marker='.')
		if not ymin:
			ymin = self.get_min(yfield, component)
		if not ymax:
			ymax = self.get_max(yfield, component)
		if not xmin:
			xmin = self.get_min(xfield)
		if not xmax:
			xmax = self.get_max(xfield)
		self.axes.set_xlim(xmin, xmax)
		self.axes.set_ylim(ymin, ymax)
		self.axes.grid(True)
		self.axes.set_title(self.title)
		# self.axes.legend()
 
	def renderS(self, xfield=field.coords,
				component = None, xmin = None, xmax = None, 
							   ymin = None, ymax = None):
		if component is not None:
			stressXX = self.storage[field.pressure]-self.storage[field.deviatorStress][:,component]
		else:
			stressXX = self.storage[field.pressure]-self.storage[field.deviatorStress][:,0]
		self.axes.plot(self.storage[xfield][:,0], stressXX , 
					label="StressXX", linestyle="", marker='.')
		if not ymin:
			ymin = min(stressXX)
		if not ymax:
			ymax = max(stressXX)
		if not xmin:
			xmin = self.get_min(xfield)
		if not xmax:
			xmax = self.get_max(xfield)
		self.axes.set_xlim(xmin, xmax)
		self.axes.set_ylim(ymin, ymax)
		self.axes.grid(True)
		self.axes.set_title(self.title)
		# self.axes.legend()

	def clear(self):
		self.axes.clear()

	def particles(self,
				 xmin = None, xmax = None, ymin = None, ymax = None, timestep = 0):
		from matplotlib import cm
		fieldList = [
		field.density,
		field.energy,
		field.pressure,
		field.energyRate,
		field.volumeRate]
		import os

		for ifield in fieldList:
			sc = self.axes.scatter(self.storage[field.coords][:,0], self.storage[field.coords][:,1], 
					c = self.storage[ifield][:],cmap=cm.jet, label=ifield)
			cb = plt.colorbar(sc)
			if not xmin:
				xmin = self.get_min(field.coords,0)
			if not xmax:
				xmax = self.get_max(field.coords,0)
			if not ymin:
				ymin = self.get_min(field.coords,1)
			if not ymax:
				ymax = self.get_max(field.coords,1)
			self.axes.set_xlim(xmin, xmax)
			self.axes.set_ylim(ymin, ymax)
			self.axes.grid(True)
			self.axes.set_title(self.title)
			dirname = "output//"+str(ifield)
			if not os.path.exists(dirname):
				os.mkdir(dirname)
			self.save(dirname+"//"+str(ifield)+"_"+str(timestep)+"_s.png", fmt='png')
			cb.remove()
			self.clear()
		fieldsPairs = [(field.force, 0),(field.force, 1), (field.velocity, 0),(field.velocity, 1),
					   (field.deviatorSRate,0)]
		for ifield, component in fieldsPairs:
			sc = self.axes.scatter(self.storage[field.coords][:,0], self.storage[field.coords][:,1], 
						c = self.storage[field.force][:,component],cmap=cm.jet, label=ifield+str(component))
			cb = plt.colorbar(sc)
			if not xmin:
				xmin = self.get_min(field.coords,0)
			if not xmax:
				xmax = self.get_max(field.coords,0)
			if not ymin:
				ymin = self.get_min(field.coords,1)
			if not ymax:
				ymax = self.get_max(field.coords,1)
			self.axes.set_xlim(xmin, xmax)
			self.axes.set_ylim(ymin, ymax)
			self.axes.grid(True)
			self.axes.set_title(self.title)
			dirname = "output//"+str(ifield)
			if not os.path.exists(dirname):
				os.mkdir(dirname)
			self.save(dirname+"//"+ifield+str(component)+"_"+str(timestep)+"_s.png", fmt='png')
			cb.remove()
			self.clear()