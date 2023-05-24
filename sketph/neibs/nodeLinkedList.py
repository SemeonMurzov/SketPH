class NodeLL:
	"""
	Nodes for a linked-list
	"""
	def __init__(self, id):
		self.id = id

		self.prev = None
		self.next = None

	def attach(self, start):
		self.removeContacts()
		# attach to the front/head
		self.next = start.next
		if self.next is not None:
			self.next.prev = self

		self.prev = start
		start.next = self
			
	def removeContacts(self):
		if self.prev is not None:
			self.prev.next = self.next
			self.prev = None
		if self.next is not None:
			self.next.prev = self.prev
			self.next = None