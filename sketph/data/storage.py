from .field import fields_list
from .field import fields_dtype
import collections.abc

import numpy as np

class Storage:    
	
	def __init__(self, fields = [], length = 0):
		self.data = dict()
		self.len = int(1.2*length)
		self.lenActive = length  
		self.fields = []
		for field in fields:
			if field in fields_list:
				self.data[field] = np.zeros(self.len, dtype=fields_dtype[field])
			else:
				raise Exception('Field ', field, ' is not found in fields list')

		for field in self.data:
			self.fields.append(field)
			
	def load(self, data):
		self.data = dict()
		self.len = 0
		self.fields = []
		k = 0 
		if type(data) == dict:
			for field in data.keys():
				if field in fields_list:
					self.data[field] = data[field]
					if self.len==0:
						self.lenActive = len(data[field])
						self.len = int(1.2*self.lenActive)
					self.fields.append(field)
				else:
					raise Exception('Field ', field, ' is not found in fields list')
			
	def add(self,nAdd):
		if self.lenActive+nAdd<=self.len:
			self.lenActive += nAdd
		else:
			self.resize(self.lenActive+nAdd,int(1.2*(self.lenActive+nAdd)))

	def delete(self,iDelete):
		lenDelete = len(iDelete)
		nDeleted = 0
		if lenDelete<self.lenActive:
			for i in iDelete:
				for field in self.fields:
					if len(self.data[field].shape)>1:
						self.data[field][i,:] = self.data[field][-(self.len-self.lenActive+1)-nDeleted,:]
					else:
						self.data[field][i] = self.data[field][-(self.len-self.lenActive+1)-nDeleted]
				nDeleted += 1
			self.lenActive -= nDeleted
		else:
			self.lenActive = 0
		return np.arange(self.lenActive,self.lenActive+nDeleted)

	def resize(self, newActiveSize, newSize):
		for field in self.fields:
			componentLen = len(self.data[field].shape)
			if componentLen>1:
				self.data[field].resize((newSize,self.data[field].shape[1]),refcheck=False)
			else:
				self.data[field].resize(newSize,refcheck=False)
		self.lenActive = newActiveSize
		self.len = newSize
			
	def as_dict(self):
		d = dict()
		for field in self.fields:
			d.setdefault(field)
			d[field] = self[field]
		return d

	def __getitem__(self, idx_field):
		if isinstance(idx_field, int) or isinstance(idx_field, np.int64):
			item = {}
			for field in self.fields:
				item[field] = self.data[field][idx_field]
			return item
		else:
			if len(np.shape(self.data[idx_field]))>1:
				return self.data[idx_field][:self.lenActive,:]
			else:
				return self.data[idx_field][:self.lenActive]
	def __setitem__(self, key, value):
		self.data[key][:self.lenActive] = value

	def __iter__ (self):
		self.n = 0
		return self
	
	def __next__(self):
		if self.n < self.len:
			item = {}
			for field in self.fields:
				item[field] = self.data[field][self.n]
			self.n += 1
			return item
		else:
			raise StopIteration