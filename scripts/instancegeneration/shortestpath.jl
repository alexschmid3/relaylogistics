
function remove!(a, item)
    deleteat!(a, findall(x->x==item, a))
end

#-----------------------------------------------------------------------#

function findshortestpath_redcost_djikstras(i, arcredcosts, numnodes, numarcs, arcs, arcLookup, A_plus, A_minus, Origin, Destination, nodesLookup, A_space, tstep)

	#Create copies of all important sets (so they can be modified for the shortest path problem)
	arcsSP, arcLookupSP, A_plusSP, nodesLookupSP = deepcopy(arcs), deepcopy(arcLookup), deepcopy(A_plus), deepcopy(nodesLookup)
	cSP = Dict()
	for a in 1:numarcs
		cSP[a] = arcredcosts[i,a]
	end

	#Create dummy origin and destination nodes + arcs 
	dummyorig = numnodes + 1
	dummydest = numnodes + 2
	numnodesSP = numnodes + 2
	A_plusSP[dummyorig], A_plusSP[dummydest] = [], []
	nodesLookupSP[dummyorig], nodesLookupSP[dummydest] = (0,0), (0,0)

	#Remove arcs to disallow leaving late from origin (after the pick up window)
	for n in Origin[i]
		n2, t2 = nodesLookup[n][1], nodesLookup[n][2] + tstep
		if t2 <= horizon
			a = arcs[n, nodes[n2, t2]]
			remove!(A_plusSP[n], a)
		end
	end

	#Disincentivize arcs that are too long to be covered by a driver
	for a in 1:numarcs
		if cSP[a] > 12
			cSP[a] = 100000
		end
	end

	#Add arcs from dummy origin to all origin nodes in pickup window
	index = numarcs + 1
	for j in 1:length(Origin[i])
		n = Origin[i][j]
		arcsSP[dummyorig, n] = index
		arcLookupSP[index] = (dummyorig, n)
		#cSP[index] = tstep * (j-1) #Cost of the arc is the time from the start of the pickup window
		cSP[index] = 0
		push!(A_plusSP[dummyorig], index)
		index += 1
	end

	#Add arcs from all destination nodes in delivery window to the dummy destination
	for n in Destination[i]
		arcsSP[n, dummydest] = index
		push!(A_plusSP[n], index)
		arcLookupSP[index] = (n, dummydest) 
		cSP[index] = 0
		index += 1
	end

	#Initialize shortest path algorithm (Dijkstra's)
	visitednodes = zeros(numnodesSP)
	currdistance = repeat([999999999.0],outer=[numnodesSP])
	currdistance[dummyorig] = 0
	currnode = dummyorig
	prevnode, prevarc = zeros(numnodesSP), zeros(numnodesSP)
	nopathexistsflag = 0
	algorithmstopflag = 0
	
	#Find shortest path from dummy origin to dummy destination
	#while (visitednodes[dummydest] == 0) & (nopathexistsflag == 0)
	while (algorithmstopflag == 0) & (nopathexistsflag == 0)

		#Assess all neighbors of current node
		for a in A_plusSP[currnode]
			n = arcLookupSP[a][2]
			if visitednodes[n] == 0

				newdist = currdistance[currnode] + cSP[a]
				
				if newdist < currdistance[n]
					currdistance[n] = newdist
					prevnode[n] = currnode
					prevarc[n] = a
				end
			end
		end

		#Mark the current node as visited
		visitednodes[currnode] = 1

		#Update the current node 
		#currdistance_unvisited = copy(currdistance)
		#for n in 1:numnodesSP
		#	if visitednodes[n] == 1
		#		currdistance_unvisited[n] = 999999999
		#	end
		#end
		currdistance_unvisited = [if visitednodes[n] == 0 currdistance[n] else 999999999 end for n in 1:numnodesSP]
		currnode = argmin(currdistance_unvisited)

		#If all remaining nodes have tentative distance of 999999999 and the algorithm has not visited the destination, then there is no path from origin to destination
		if (minimum(currdistance_unvisited) == 999999999) & (visitednodes[dummydest] == 0)
			nopathexistsflag = 1
		elseif minimum(currdistance_unvisited) == 999999999
			algorithmstopflag = 1
		end

	end

	#if nopathexistsflag == 1
	#	println("==================================================")
	#	println("==================================================")
	#	println("nopathexistsflag = ", nopathexistsflag, " for order $i")
	#	println("==================================================")
	#	println("==================================================")
	#end

	#Format the shortest path output
	patharcs = [0 for a in 1:numarcs]
	node = Int(prevnode[dummydest])

	while Int(prevnode[node]) != dummyorig
		patharcs[Int(prevarc[node])] = 1
		node = Int(prevnode[node])
	end

	#Format the shortest path output
	shortestpathnodes_rev = [dummydest]
	shortestpatharcs_rev = []
	node = dummydest
	while node != dummyorig
		push!(shortestpatharcs_rev, Int(prevarc[node]))
		node = Int(prevnode[node])
		push!(shortestpathnodes_rev, node)
	end
	shortestpathnodes = reverse(shortestpathnodes_rev[2:length(shortestpathnodes_rev)-1]) 
	shortestpatharcs = reverse(shortestpatharcs_rev[2:length(shortestpatharcs_rev)-1]) 

	return currdistance[dummydest], shortestpathnodes, shortestpatharcs

end

#-----------------------------------------------------------------------#

function findshortestpath_redcost_multiples(i, p, arcredcosts, numnodes, numarcs, extendednumarcs, arcs, arcLookup, Origin, Destination, nodesLookup, tstep)

	#i, p, arcredcosts, numnodes, numarcs, arcs, arcLookup, Origin, Destination, nodesLookup, tstep = i, pathsperiter, arcredcosts, extendednumnodes, extendednumarcs, extendedarcs, arcLookup, Origin, Destination, nodesLookup, tstep
	
	#Create copies of all important sets (so they can be modified for the shortest path problem)
	arcsSP, arcLookupSP, nodesLookupSP = deepcopy(arcs), deepcopy(arcLookup), deepcopy(nodesLookup)
	cSP = Dict()
	for a in 1:extendednumarcs
		cSP[a] = arcredcosts[i,a]
	end

	#Create dummy origin and destination nodes + arcs 
	dummyorig = numnodes + 1
	dummydest = numnodes + 2
	numnodesSP = numnodes + 2
	nodesLookupSP[dummyorig], nodesLookupSP[dummydest] = (0,0), (0,0)

	#Add arcs from dummy origin to all origin nodes in pickup window
	index = extendednumarcs + 1
	for j in 1:length(Origin[i])
		n = Origin[i][j]
		arcsSP[dummyorig, n] = index
		arcLookupSP[index] = (dummyorig, n)
		#cSP[index] = tstep * (j-1) #Cost of the arc is the time from the start of the pickup window
		cSP[index] = 0
		index += 1
	end

	#Add arcs from all destination nodes in delivery window to the dummy destination
	for n in Destination[i]
		arcsSP[n, dummydest] = index
		arcLookupSP[index] = (n, dummydest) 
		cSP[index] = 0
		index += 1
	end
	numarcsSP = length(cSP)
	arclistSP = [a for a in 1:length(cSP)]

	#Remove arcs to disallow leaving late from origin (after the pick up window)
	if !(i in ordersinprogress)
		for n in Origin[i]
			n2, t2 = nodesLookup[n][1], nodesLookup[n][2] + tstep
			if t2 <= horizon
				a = arcs[n, nodes[n2, t2]]
				remove!(arclistSP, a)
				#println(nodesLookup[arcLookup[a][1]] , " ==> ", nodesLookup[arcLookup[a][2]] )
			end
		end
	end

	#Remove arcs to waiting around at the destination (during the delivery window)
	for n in Destination[i]
		n2, t2 = nodesLookup[n][1], nodesLookup[n][2] + tstep
		if t2 <= horizon
			a = arcs[n, nodes[n2, t2]]
			remove!(arclistSP, a)
			#println(nodesLookup[arcLookup[a][1]] , " ==> ", nodesLookup[arcLookup[a][2]] )
		end
	end

	#Remove arcs arriving at origin 
	for n in Origin[i], a in A_minus[n]
		remove!(arclistSP, a)
		#println(nodesLookup[arcLookup[a][1]] , " ==> ", nodesLookup[arcLookup[a][2]] )
	end

	#Remove arcs leaving from destination 
	for n in Destination[i], a in A_plus[n]
		remove!(arclistSP, a)
		#println(nodesLookup[arcLookup[a][1]] , " ==> ", nodesLookup[arcLookup[a][2]] )
	end

	#Remove arcs that are not feasible because they are longer than the driver shift
	for a in 1:numarcs
		arctime = nodesLookup[arcLookup[a][2]][2] - nodesLookup[arcLookup[a][1]][2]
		if arctime > shiftlength
			remove!(arclistSP, a)
		end
	end

	#Initialize shortest path algorithm (Bellman-Ford)
	currdistance = repeat([999999999.0],outer=[numnodesSP])
	currdistance[dummyorig] = 0
	prevnode, prevarc = zeros(numnodesSP), zeros(numnodesSP)
	orderstarttime = nodesLookup[Origin[i][1]][2]
	maxspiter = (horizon - orderstarttime)/tstep + 3
	
	for iteration in 1:maxspiter #Max path length
		for a in arclistSP #1:numarcsSP
			n_end, n_start = arcLookupSP[a][2], arcLookupSP[a][1]
			if currdistance[n_end] > currdistance[n_start] + cSP[a] + 0.000001
				currdistance[n_end] = currdistance[n_start] + cSP[a]
				prevnode[n_end] = n_start
				prevarc[n_end] = a
			end
		end
	end

	#Format the shortest path output
	shortestpathnodes_rev = [dummydest]
	shortestpatharcs_rev = []
	node = dummydest
	while node != dummyorig
		push!(shortestpatharcs_rev, Int(prevarc[node]))
		node = Int(prevnode[node])
		push!(shortestpathnodes_rev, node)
	end
	shortestpathnodes = reverse(shortestpathnodes_rev[2:length(shortestpathnodes_rev)-1]) 
	shortestpatharcs = reverse(shortestpatharcs_rev[2:length(shortestpatharcs_rev)-1]) 

	#Find up to p-1 additional paths
	numpaths = 1
	secondarypathnodes, secondarypatharcs, secondarypathredcost = [], [], []
	shortestdestnode = prevnode[dummydest]
	if p < numpaths
		for n in setdiff(Destination[i], shortestdestnode)
			if currdistance[n] < -0.000001
				node = n
				newpath = []
				while node != dummyorig
					push!(secondarypathnodes, node)
					if prevnode[node] != dummyorig
						push!(newpath, Int(prevarc[node]))
					end
					node = Int(prevnode[node])
				end
				push!(secondarypatharcs, reverse(newpath))
				push!(secondarypathredcost, currdistance[n])
				numpaths += 1
				if numpaths >= p
					break
				end
			end
		end
	end

	return currdistance[dummydest], shortestpathnodes, shortestpatharcs, numpaths, secondarypathnodes, secondarypatharcs, secondarypathredcost

end

#-----------------------------------------------------------------------#

function findshortestpath_redcost_bellmanford(i, arcredcosts, numnodes, numarcs, arcs, arcLookup, Origin, Destination, nodesLookup, tstep)

	#Create copies of all important sets (so they can be modified for the shortest path problem)
	arcsSP, arcLookupSP, nodesLookupSP = deepcopy(arcs), deepcopy(arcLookup), deepcopy(nodesLookup)
	cSP = Dict()
	for a in 1:numarcs
		cSP[a] = arcredcosts[i,a]
	end

	#Create dummy origin and destination nodes + arcs 
	dummyorig = numnodes + 1
	dummydest = numnodes + 2
	numnodesSP = numnodes + 2
	nodesLookupSP[dummyorig], nodesLookupSP[dummydest] = (0,0), (0,0)

	#Add arcs from dummy origin to all origin nodes in pickup window
	index = numarcs + 1
	for j in 1:length(Origin[i])
		n = Origin[i][j]
		arcsSP[dummyorig, n] = index
		arcLookupSP[index] = (dummyorig, n)
		cSP[index] = tstep * (j-1) #Cost of the arc is the time from the start of the pickup window
		index += 1
	end

	#Add arcs from all destination nodes in delivery window to the dummy destination
	for n in Destination[i]
		arcsSP[n, dummydest] = index
		arcLookupSP[index] = (n, dummydest) 
		cSP[index] = 0
		index += 1
	end

	numarcsSP = length(cSP)
	arclistSP = [a for a in 1:length(cSP)]

	#Remove arcs to disallow leaving late from origin (after the pick up window)
	for n in Origin[i]
		n2, t2 = nodesLookup[n][1], nodesLookup[n][2] + tstep
		if t2 <= horizon
			a = arcs[n, nodes[n2, t2]]
			remove!(arclistSP, a)
		end
	end

	#Remove arcs to waiting around at the destination (during the delivery window)
	for n in Destination[i]
		n2, t2 = nodesLookup[n][1], nodesLookup[n][2] + tstep
		if t2 <= horizon
			a = arcs[n, nodes[n2, t2]]
			remove!(arclistSP, a)
		end
	end

	#Remove arcs that are not feasible because they are longer than the driver shift
	for a in 1:numarcs
		arctime = nodesLookup[arcLookup[a][2]][2] - nodesLookup[arcLookup[a][1]][2]
		if arctime > shiftlength
			remove!(arclistSP, a)
		end
	end

	#Initialize shortest path algorithm (Bellman-Ford)
	currdistance = repeat([999999999.0],outer=[numnodesSP])
	currdistance[dummyorig] = 0
	prevnode, prevarc = zeros(numnodesSP), zeros(numnodesSP)
	orderstarttime = nodesLookup[Origin[i][1]][2]
	maxspiter = (horizon - orderstarttime)/tstep + 3
	
	for iteration in 1:maxspiter #Max path length
		for a in arclistSP 
			n_end, n_start = arcLookupSP[a][2], arcLookupSP[a][1]
			if currdistance[n_end] > currdistance[n_start] + cSP[a] + .000001
				currdistance[n_end] = currdistance[n_start] + cSP[a]
				prevnode[n_end] = n_start
				prevarc[n_end] = a
			end
		end
	end

	#Format the shortest path output
	patharcs = [0 for a in 1:numarcs]
	node = Int(prevnode[dummydest])

	while Int(prevnode[node]) != dummyorig
		patharcs[Int(prevarc[node])] = 1
		node = Int(prevnode[node])
	end

	#Format the shortest path output
	shortestpathnodes_rev = [dummydest]
	shortestpatharcs_rev = []
	node = dummydest
	while node != dummyorig
		push!(shortestpatharcs_rev, Int(prevarc[node]))
		node = Int(prevnode[node])
		push!(shortestpathnodes_rev, node)
	end
	shortestpathnodes = reverse(shortestpathnodes_rev[2:length(shortestpathnodes_rev)-1]) 
	shortestpatharcs = reverse(shortestpatharcs_rev[2:length(shortestpatharcs_rev)-1]) 
	
	return currdistance[dummydest], patharcs, shortestpathnodes, shortestpatharcs

end

#-----------------------------------------------------------------------#

function findshortestpath_redcost_pbcg(i, arcredcosts, numnodes, numarcs, extendednumarcs, arcs, nodes)
										
	#Create copies of all important sets (so they can be modified for the shortest path problem)
	arcsSP, arcLookupSP, nodesLookupSP = deepcopy(arcs), deepcopy(arcLookup), deepcopy(nodesLookup)
	cSP = Dict()
	for a in 1:extendednumarcs
		cSP[a] = arcredcosts[i,a]
	end

	#Create dummy origin and destination nodes + arcs 
	dummyorig = numnodes + 1
	dummydest = numnodes + 2
	numnodesSP = numnodes + 2
	nodesLookupSP[dummyorig], nodesLookupSP[dummydest] = (0,0), (0,0)

	#Add arcs from dummy origin to all origin nodes in pickup window
	index = extendednumarcs + 1
	for j in 1:length(Origin[i])
		n = Origin[i][j]
		arcsSP[dummyorig, n] = index
		arcLookupSP[index] = (dummyorig, n)
		#cSP[index] = tstep * (j-1) 		#Cost of the arc is the time from the start of the pickup window
		cSP[index] = 0
		index += 1
	end

	#Add arcs from all destination nodes in delivery window to the dummy destination
	for n in Destination[i]
		arcsSP[n, dummydest] = index
		arcLookupSP[index] = (n, dummydest) 
		cSP[index] = 0
		index += 1
	end

	#Create arc list
	numarcsSP = length(cSP)
	arclistSP = [a for a in 1:length(cSP)]

	#Remove arcs to disallow leaving late from origin (after the pick up window)
	if !(i in ordersinprogress)
		for n in Origin[i]
			n2, t2 = nodesLookup[n][1], nodesLookup[n][2] + tstep
			if t2 <= horizon
				a = arcs[n, nodes[n2, t2]]
				remove!(arclistSP, a)
			end
		end
	end

	#Remove arcs to waiting around at the destination (during the delivery window)
	for n in Destination[i]
		n2, t2 = nodesLookup[n][1], nodesLookup[n][2] + tstep
		if t2 <= horizon
			a = arcs[n, nodes[n2, t2]]
			remove!(arclistSP, a)
		end
	end

	#Remove arcs arriving at origin 
	for n in Origin[i], a in A_minus[n]
		remove!(arclistSP, a)
	end

	#Remove arcs leaving from destination 
	for n in Destination[i], a in A_plus[n]
		remove!(arclistSP, a)
	end

	#Remove arcs that are not feasible because they are longer than the driver shift
	for a in 1:numarcs
		arctime = nodesLookup[arcLookup[a][2]][2] - nodesLookup[arcLookup[a][1]][2]
		if arctime > shiftlength
			remove!(arclistSP, a)
		end
	end

	#Initialize shortest path algorithm (Bellman-Ford)
	currdistance = repeat([999999999.0],outer=[numnodesSP])
	currdistance[dummyorig] = 0
	prevnode, prevarc = zeros(numnodesSP), zeros(numnodesSP)
	orderstarttime = nodesLookup[Origin[i][1]][2]
	maxspiter = (horizon - orderstarttime)/tstep + 3
	
	for iteration in 1:maxspiter #Max path length
		for a in arclistSP 
			n_end, n_start = arcLookupSP[a][2], arcLookupSP[a][1]
			if currdistance[n_end] > currdistance[n_start] + cSP[a] + .000001
				currdistance[n_end] = currdistance[n_start] + cSP[a]
				prevnode[n_end] = n_start
				prevarc[n_end] = a
			end
		end
	end

	#Format the shortest path output
	patharcs = [0 for a in 1:extendednumarcs]
	node = Int(prevnode[dummydest])

	while Int(prevnode[node]) != dummyorig
		patharcs[Int(prevarc[node])] = 1
		node = Int(prevnode[node])
	end

	#Format the shortest path output
	shortestpathnodes_rev = [dummydest]
	shortestpatharcs_rev = []
	node = dummydest
	while node != dummyorig
		push!(shortestpatharcs_rev, Int(prevarc[node]))
		node = Int(prevnode[node])
		push!(shortestpathnodes_rev, node)
	end
	shortestpathnodes = reverse(shortestpathnodes_rev[2:length(shortestpathnodes_rev)-1]) 
	shortestpatharcs = reverse(shortestpatharcs_rev[2:length(shortestpatharcs_rev)-1]) 
	
	return currdistance[dummydest], patharcs, shortestpathnodes, shortestpatharcs

end

#-----------------------------------------------------------------------#

