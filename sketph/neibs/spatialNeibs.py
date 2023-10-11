import numpy as np
from .nodeLinkedList import NodeLL
from sketph.box.rectangle import Box
from sketph.data import field
from sketph.data import Storage
dim = 2
import itertools

class SpatialNeibs:
	"""
	class for neibs detection in the computational box
	we use a hash-table which is consist of cells with linked list, where each particle is a NodeLL
	the lists rebuild each time the size of the mesh is changed
	"""
	def __init__(self,
			box = Box(),
			particles  = Storage(fields = [field.coords, field.size]), 
			smoothingScale = 1.0,
			buf = 2.8):
		self.IndexAxis = {0:"x",1:"y"}
		self.box = box
		self.ss  = smoothingScale
		# initialize the cells and fill it with particles 
		# calculate the minimal cell size to cover the box where particles stored 
		self.buf = buf
		self.refCellSize = buf*np.max(particles[field.size]*self.ss)
		self.adjustCells()
		# create empty cells of linked lists of nodes - particles indexes
		self.nodes = {}
		self.createCells()
		self.fillCells(particles)
		self.setLocalCellsIndexes()
		print("spatial neibs table")
		print(self.allCellsNumber, "all Cells number ") 
		print(self.cellsNumber[0], "all Cells number X") 
		print(self.cellsNumber[1], "all Cells number Y") 

	# set dictionary of nodes empty and initialize cells
	def createCells(self):
		self.celListNodes = [NodeLL(-1) for _ in range(self.allCellsNumber)]

	# set new celListNodes
	def reCreate(self):
		del self.celListNodes
		self.createCells()

	# fill cells NodeLL objects linking them with specified cells linked list
	def fillCells(self,particles):
		for i in self.celListNodes:
			i.removeContacts()
		for i in range(particles.lenActive):
			self.moveNode(i, particles[field.coords][i])

	# attach node associated with id to the Cell LL anyway	
	def moveNode(self, id, new_coords):
		iCell = self.whichCell(new_coords)
		if not id in self.nodes:
	# create a new NodeLL for particle and  add to the global dictionary of nodes
			node = NodeLL(id)
			self.nodes[id] = node
			node.attach(self.celListNodes[iCell])
		else:
			if self.nodes[id].id!=id:
				self.nodes[id].id = id
			self.nodes[id].attach(self.celListNodes[iCell])

	# replace node linked List Contacts to erase it from the SPH sum	
	def replaceNodes(self,idsRemove,idsReplace):
						#  ideleted,idelete[0]
		for i in range(len(idsRemove)):
			self.nodes[idsReplace[i]].removeContacts()
			swap = self.nodes[idsReplace[i]]
			self.nodes[idsReplace[i]] = self.nodes[idsRemove[i]]
			self.nodes[idsReplace[i]].id = idsReplace[i]
			self.nodes[idsRemove[i]] = swap

	# calculate index of cell where the particle is situated
	def getIndexSlice(self,k):
		return 1 if k==0 else self.cellsNumber[k-1]*self.getIndexSlice(k-1)
	def whichCell(self,coords):
		iCell=0
		for k in range(dim):
			iCell += int((coords[k]-self.box.Lmin[k])/self.actCellSize[k])*self.getIndexSlice(k)
		return iCell
	""" adjusting the size of cells to the box size for a given current minimal size of cell
		We link particles through the periodic boundary to get neighbours using real periodic boundary
		and avoid this quantization inconsistency of the hash-table and dimensions of the simulation box.
	"""
	def adjustCells(self):
		self.cellsNumber  = np.zeros(dim,dtype='i4')
		for k in range(dim):
			self.cellsNumber[k] = int(self.box.L[k]/self.refCellSize)
		self.allCellsNumber  = np.prod(self.cellsNumber)
		self.actCellSize = np.zeros(dim, dtype='f8')
		for k in range(dim):
			self.actCellSize[k] = self.box.L[k]/self.cellsNumber[k]
	""" set a local neighbours indexes of the cells in the rectangular
		mesh with a periodic boundary conditions for the internal cells
	"""
	def period(self,iDirectCell,k):
		if iDirectCell>self.cellsNumber[k]-1:
			return iDirectCell-self.cellsNumber[k]
		if iDirectCell<0:
			return iDirectCell+self.cellsNumber[k]
		return iDirectCell
	# for _, _ in itertools.product(range(3), range(3)): gives that order
	# 0 0, 0 1, 0 2, 1 0, 1 1, 1 2, 2 0, 2 1, 2 2
	# set all local cells neighbours including periodic boundary conditions
	def setLocalCellsIndexes(self):
		self.CellsNeibsIndexes = [[[-1] for _ in range(3*3)]
										for _ in range(self.cellsNumber[0]*self.cellsNumber[1])]
		for iXCell, iYCell in itertools.product(range(self.cellsNumber[0]),
												range(self.cellsNumber[1])):
			iCell = iXCell+iYCell*self.cellsNumber[0]
			for dX, dY in itertools.product(range(-1,2), range(-1,2)):
				# periodic boundary conditions including side case and corner case
				iXNeib = self.period(iXCell+dX,0)
				iYNeib = self.period(iYCell+dY,1)
				self.CellsNeibsIndexes[iCell][dX+1+3*(dY+1)] = iXNeib+iYNeib*self.cellsNumber[0]

	""" make a full rebuild of cells, because the size of the mesh is too small
		if the size of particles is consistent, cells are refilled by
		the current configuration of particles
	"""
	def checkUpdate(self, particles, buf = 1.0):
		expectCellSize = buf*np.max(particles[field.size]*self.ss)
		# build new cells 
		if self.refCellSize<expectCellSize:
			self.refCellSize = expectCellSize
			self.adjustCells()
			self.createCells()
			self.fillCells(particles)
			self.setLocalCellsIndexes()
		else:
			self.fillCells(particles)
	""" get all neibs through generator for given
		point with coords[x],coords[y] 
	"""
	def getIt(self, coords):
		iCell = self.whichCell(coords)
		for loc in range(len(self.CellsNeibsIndexes[iCell][:])):
			ineibCell = self.CellsNeibsIndexes[iCell][loc]
			node = self.celListNodes[ineibCell].next
			while node is not None:
				yield node.id
				node = node.next


