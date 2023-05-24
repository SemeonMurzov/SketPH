from sketph.data import field
import numpy as np
def periodic_indexing(start,stop,length):
	if start <= stop:
		return (np.arange(start,stop),)
	else:
		return (np.concatenate([np.arange(start,length),np.arange(0,stop)]),)

def feed(particles, neibsTable, upIn, xInsert, inflowSample = None, 
		  firstInsertIndex = None, xStart = None, 
		  boxInflowPeriod = None, inflowWindowPosition = None):
	# index of first particle to be inserted
	iStart = firstInsertIndex
	# distance from right boundary x = xInsert to the plane x = inflowWindowPosition
	lengthToFill = inflowWindowPosition-xInsert
	# periodic indexes
	Lx = boxInflowPeriod
	inflowSample[field.coords][0:iStart,0] += Lx
	''' index of the first particle to be inserted on the next interation of AMW 
	    which can be lesser than current first index
	 	find iEnd: inflowSample[iEnd-1]<=xStart+lengthToFill<inflowSample[iEnd] '''
	iEnd = np.searchsorted(inflowSample[field.coords][iStart:,0],xStart+lengthToFill,side='right')+iStart
	xStartNewPeriod = 0
	# correcting indexes due to the periodicity of feeding
	if iEnd >= inflowSample.lenActive:
		iEnd = np.searchsorted(inflowSample[field.coords][:iStart,0],xStart+lengthToFill,side='right')
		xStartNewPeriod = inflowSample[field.coords][iEnd-1,0]-Lx
	# form a numpy array range of indexes
	inflowParticlesIndexes = periodic_indexing(iStart,iEnd,inflowSample.lenActive)
	# change the size of the storage to fit into the new number of particles 
	nParticlesToInsert = len(inflowParticlesIndexes[0])
	if nParticlesToInsert>0:
		oldSize = particles.lenActive
		if oldSize + nParticlesToInsert>particles.len:
			particles.resize(oldSize + nParticlesToInsert, int(1.2*(oldSize + nParticlesToInsert)))
		else:
			particles.lenActive += nParticlesToInsert
		iInsert = np.arange(oldSize,oldSize + nParticlesToInsert)
		# insert particles to the storage accprding to the coords shift

		particles[field.coords][iInsert,0] = \
					inflowSample[field.coords][inflowParticlesIndexes,0]-xStart+xInsert
		# boundary of the sample is a coordinate of the last inserted particle
		xInsert = particles[field.coords][-1,0]
		# go back within periodic boundary conditions for the inserting sample
		inflowSample[field.coords][0:iStart,0] -= Lx
		# track the position of the last inserted particle in it's local coordinates
		if xStartNewPeriod:
			xStart = xStartNewPeriod
		else:
			xStart = inflowSample[field.coords][iEnd-1,0]
		particles[field.coords][iInsert,1] = inflowSample[field.coords][inflowParticlesIndexes,1]
		particles[field.material][iInsert] = inflowSample[field.material][inflowParticlesIndexes]
		particles[field.density][iInsert]  = inflowSample[field.density][inflowParticlesIndexes]
		particles[field.mass][iInsert]     = inflowSample[field.mass][inflowParticlesIndexes]
		particles[field.energy][iInsert]   = inflowSample[field.energy][inflowParticlesIndexes]
		particles[field.size][iInsert]     = inflowSample[field.size][inflowParticlesIndexes]
		particles[field.velocity][iInsert,0] = upIn
		particles[field.velocity][iInsert,1] = 0    
		particles[field.deviatorStress][iInsert,:] = 0    
		particles[field.pressure][iInsert] = 0
		# create new nodes for neibsTable
		for id in iInsert:
			neibsTable.moveNode(id, particles[field.coords][id])
	
		print("nParticleToInsert, new insert position ", nParticlesToInsert,particles[field.coords][iInsert[-1],0])
		return xStart, iEnd, xInsert
	else:
		# return the coords after periodicity
		inflowSample[field.coords][0:iStart,0] -= Lx
		return xStart, firstInsertIndex, xInsert