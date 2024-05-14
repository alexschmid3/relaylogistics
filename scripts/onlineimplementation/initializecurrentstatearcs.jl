
function finddriversets_online(T_off_online, driverStartNodes, lasttimehome)
		
	#Find the effective shifts 
	#i.e. For a three day horizon starting on a Monday, a MTWTF driver looks the same as a MTWR+Su driver
	effectiveshifts = unique!(T_off_online)
	numeffshifts = length(effectiveshifts)
	effshift, shiftsincluded = [], [[] for es in 1:numeffshifts]
	for s in 1:numshifts
		es = findfirst(x->x==T_off_online[s], effectiveshifts)
		push!(effshift, es)
		push!(shiftsincluded[es], s)
	end

	#Create driver sets
	driversingroup = Dict()
	driversets = []
	for homeloc in 1:numlocs, shiftsched in 1:numshifts, lth in -1*10*24:tstep:10*24, startnode in 1:numnodes
		driversingroup[homeloc, shiftsched, startnode, lth] = []
	end
	for d in drivers
		push!(driversingroup[driverHomeLocs[d], drivershift[d], driverStartNodes[d], lasttimehome[d]], d)
	end
	for item in driversingroup
		if item[2] != []
			push!(driversets, item[1])
		end
	end

	#Group lookups
	drivergroupnum, drivergroupdesc = Dict(), Dict()
	index = 1
	for (hl,ss,sn,lth) in driversets
		drivergroupnum[hl,ss,sn,lth] = index
		drivergroupdesc[index] = (hl,ss,sn,lth)
		index += 1
	end
	numdrivergroups = length(drivergroupdesc)

	return driversets, driversingroup, numdrivergroups, drivergroupnum, drivergroupdesc, numeffshifts, effshift, shiftsincluded

end

#----------------------------------------------------------------------------------------#

function addnextsegment_online(g, ss, currentfrags, completedfrags, deadline, offtimes, homeloc, driversingroup, drivergroupdesc, driverarcs)
	
	currfragnum = length(currentfrags)

	#println(keys(driverarcs))

	#Note this does not work for irregular work schedules!!
	offtimes_ext = union(offtimes, [t+24 for t in intersect([t for t in horizon-3*tstep:tstep:horizon],offtimes)])
	T_on_0_ext = [t for t in setdiff(0:tstep:horizon+24, offtimes_ext) if !(t-tstep in setdiff(0:tstep:horizon+24, offtimes_ext))]
	offtimes_0 = [t for t in offtimes_ext if !(t-tstep in offtimes_ext)]	

	#Add next segment to each fragment
	for f in 1:currfragnum
		frag = popfirst!(currentfrags)
		currloc, currtime = last(frag)
		d_ex = driversingroup[drivergroupdesc[g]][1]
		endofcurrshift = minimum(union([t for t in offtimes_0 if t > currtime], horizon+24))
		startofnextshift = minimum(union([t for t in T_on_0_ext if t > currtime], horizon+24))
		#println("   $currloc, $currtime")
		for a in driverarcs.A_plus[d_ex, nodes[currloc, currtime]]
			arcendloc, arcendtime = nodesLookup[arcLookup[a][2]]
			if traveltimebetweenlocs_rdd[arcendloc, homeloc] > 0
				effectivedeadline = maximum([t2 for t2 in offtimes_0 if t2 <= deadline])
			else
				effectivedeadline = deadline
			end
			if (arcendtime <= min(effectivedeadline - traveltimebetweenlocs_rdd[arcendloc, homeloc], endofcurrshift)) & (!((arcendloc == homeloc) & (currtime in offtimes_ext) & (setdiff([t for t in startofnextshift:tstep:deadline-tstep], offtimes) != [])) || startofnextshift > horizon)
				newfrag = deepcopy(frag)
				push!(newfrag, (arcendloc, arcendtime))
				if (arcendtime == deadline) || (arcendtime == horizon)
					push!(completedfrags, newfrag)
				else 
					push!(currentfrags, newfrag)
				end
			end
		end
	end

	return currentfrags, completedfrags

end

#----------------------------------------------------------------------------------------#

function createghostTSN(maxnightsaway)
	
	newhorizon = horizon + 2*24*maxnightsaway

	#Nodes
	newnodes = copy(nodes)
	newnodesLookup = copy(nodesLookup)
	index = numnodes + 1
	for t in horizon+tstep:tstep:newhorizon
		for l in 1:numlocs
			newnodes[l,t] = index
			newnodesLookup[index] = (l,t)
			index += 1
		end	
	end
	newnumnodes = length(newnodes)

	#Arcs
	newarcs, newarcLookup, newA_plus, newA_minus = Dict(), Dict(), Dict(), Dict()
	for node in 1:newnumnodes
		newA_plus[node] = []
		newA_minus[node] = []
	end
	index = 1
	for arc in prearcs
		traveltime = arc[3]
		if arc[1] != arc[2]
			for t in 0:tstep:newhorizon - traveltime
				startnode = newnodes[arc[1],t]
				endnode = newnodes[arc[2],t+arc[3]]
				newarcs[(startnode,endnode)] = index
				newarcLookup[index] = (startnode,endnode)	
				push!(newA_plus[startnode], index)
				push!(newA_minus[endnode], index)
				index += 1
			end
		end
	end
	for loc in 1:numlocs, t in 0:tstep:newhorizon - tstep
		startnode = newnodes[loc,t]
		endnode = newnodes[loc,t + tstep]
		newarcs[(startnode,endnode)] = index
		newarcLookup[index] = (startnode,endnode)
		push!(newA_plus[startnode], index)
		push!(newA_minus[endnode], index)	
		index += 1
	end
	newnumarcs = length(newarcs)

	#Arcs between (for shortestpathonTSN function)
	newarcsbetween = Dict()
	sortedarclist = sort([a for a in 1:newnumarcs], by=x->newnodesLookup[newarcLookup[x][1]][2])
	for t1 in 0:tstep:newhorizon, t2 in t1:tstep:newhorizon
		newarcsbetween[t1,t2] = [a for a in sortedarclist if (newnodesLookup[newarcLookup[a][1]][2] >= t1) & (newnodesLookup[newarcLookup[a][2]][2] <= t2)]
	end
	newarcsbetween_back = Dict()
	sortedarclist_back = reverse(sort([a for a in 1:newnumarcs], by=x->newnodesLookup[newarcLookup[x][2]][2]))
	for t1 in 0:tstep:newhorizon, t2 in t1:tstep:newhorizon
		newarcsbetween_back[t1,t2] = [a for a in sortedarclist_back if (newnodesLookup[newarcLookup[a][1]][2] >= t1) & (newnodesLookup[newarcLookup[a][2]][2] <= t2)]
	end

	#for i in 1:100:length(myarcs)
	#	println(newnodesLookup[newarcLookup[myarcs[i]][2]][2])
	#end

	ghosttsn = (arcsbetween=newarcsbetween, arcsbetween_back=newarcsbetween_back, numlocs=numlocs, 
			arcLookup=newarcLookup, nodesLookup=newnodesLookup, nodes=newnodes, numarcs=newnumarcs, 
			numnodes=newnumnodes, horizon=newhorizon, tstep=tstep, arcs=newarcs, extendednumarcs=newnumarcs, 
			extendednumnodes=newnumnodes, A_minus=newA_minus, A_plus=newA_plus)

	#=include("scripts/onlineimplementation/initializecurrentstatearcs.jl")
	include("scripts/visualizations/timespacenetwork.jl")
	myarcs = intersect(tsn.arcsbetween[t1,t2], setofallowedarcs)
	timespacenetwork_providetsn(ghosttsn, "outputs/viz/aaaaghost.png", [truearcs], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
	=#
	return ghosttsn

end

#----------------------------------------------------------------------------------------#

function findshortestjourneyhome(startnode, hl, deadlinetime, setofallowedarcs, tsn)
	
	deadlinenode = tsn.nodes[hl, deadlinetime]
	t1,t2 = tsn.nodesLookup[startnode][2], tsn.nodesLookup[deadlinenode][2]
	
	#Initialize shortest path algorithm 
	currdistance = repeat([999999999.0],outer=[tsn.numnodes])
	currdistance[startnode] = 0
	#currarcdistance = repeat([999999999.0],outer=[tsn.numarcs])
	prevarc = repeat([-1],outer=[tsn.numnodes])
	
	#Loop over time-space arcs in order of start time
	for a in intersect(tsn.arcsbetween[t1,t2], setofallowedarcs)
		n_start, n_end = tsn.arcLookup[a]
		if currdistance[n_start] <= 1e-4 
			currdistance[n_end] = currdistance[n_start] 
			prevarc[n_end] = a
		end
	end

	#Find the path
	for t in t1:tstep:t2
		n_end = tsn.nodes[hl, t]
		if currdistance[n_end] < 1e-4
			#Found shortest path back home
			patharcs = []
			currnode = n_end
			while currnode != startnode
				a = prevarc[currnode]
				push!(patharcs, a)
				currnode = tsn.arcLookup[a][1]
			end
			return patharcs, t
		end	
	end

	return [], dummyendtime

end

#----------------------------------------------------------------------------------------#

function shortestpathonTSN(startnode, hl, deadlinetime, setofallowedarcs, tsn)
	
	deadlinenode = tsn.nodes[hl, deadlinetime]
	t1,t2 = tsn.nodesLookup[startnode][2], tsn.nodesLookup[deadlinenode][2]
	
	#Initialize shortest path algorithm 
	currdistance = repeat([999999999.0],outer=[tsn.numnodes])
	currdistance[startnode] = 0
	currarcdistance = repeat([999999999.0],outer=[tsn.numarcs])
	
	#Loop over time-space arcs in order of start time
	for a in intersect(tsn.arcsbetween[t1,t2], setofallowedarcs)
		n_start, n_end = tsn.arcLookup[a]
		if currdistance[n_start] <= 1e-4 
			currdistance[n_end] = currdistance[n_start] 
			currarcdistance[a] = currdistance[n_start]
		end
	end

	#Run through backwards
	currdistance_back = repeat([999999999.0],outer=[tsn.numnodes])
	currdistance_back[deadlinenode] = 0
	currarcdistance_back = repeat([999999999.0],outer=[tsn.numarcs])
	for a in intersect(tsn.arcsbetween_back[t1,t2], setofallowedarcs)
		n_start, n_end = tsn.arcLookup[a]
		if currdistance_back[n_end] <= 1e-4
			currdistance_back[n_start] = currdistance_back[n_end] 
			currarcdistance_back[a] = currdistance[n_end]
		end
	end

	frontarcs = [a for a in 1:tsn.numarcs if currarcdistance[a] < 1e-4]
	backarcs = [a for a in 1:tsn.numarcs if currarcdistance_back[a] < 1e-4]
	truearcs = intersect(frontarcs, backarcs)
	#timespacenetwork("outputs/viz/aaa_all.png", [truearcs], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)

	return truearcs

end

#----------------------------------------------------------------------------------------#

function formatarcset(truearcs, tsn)

	A_plus_j = Dict()
	for n in 1:tsn.numnodes
		A_plus_j[n] = []
	end
	for a in truearcs
		n = tsn.arcLookup[a][1]
		push!(A_plus_j[n], a)
	end

	journeyarcs = (A=truearcs, A_plus=A_plus_j)

	return journeyarcs

end

#----------------------------------------------------------------------------------------#

function enumeratealljourneys(sn, journeyendingnodelist, journeyarcs)

	currentfragments = [[(sn,0)]]
	completedjourneys = []

	while currentfragments != []
		frag = pop!(currentfragments)
		for a in journeyarcs.A_plus[last(frag)[1]]
			newfrag = copy(frag)
			n = arcLookup[a][2]
			push!(newfrag, (n,a))
			if n in journeyendingnodelist
				push!(completedjourneys, [item[2] for item in newfrag[2:length(newfrag)]])
			else
				push!(currentfragments, newfrag)
			end
		end
	end

	return completedjourneys

end

#----------------------------------------------------------------------------------------# 

function translatearcsbetweenTSNs(arcset, originaltsn, targettsn)

	newarcset = []

	for a in arcset
		l1, t1 = originaltsn.nodesLookup[originaltsn.arcLookup[a][1]]
		l2, t2 = originaltsn.nodesLookup[originaltsn.arcLookup[a][2]]
		a_prime = targettsn.arcs[targettsn.nodes[l1,t1], targettsn.nodes[l2,t2]]
		push!(newarcset, a_prime)
	end

	return newarcset 

end

#----------------------------------------------------------------------------------------# 

function translatearcsbetweenTSNs_withextendedarcs(arcset, originaltsn, targettsn)

	newarcset = []

	for a in arcset
		l1, t1 = originaltsn.nodesLookup[originaltsn.arcLookup[a][1]]
		l2, t2 = originaltsn.nodesLookup[originaltsn.arcLookup[a][2]]
		if t1 >= targettsn.horizon     #Remove arcs that begin on or after horizon
			1+1
		elseif t2 > targettsn.horizon  #Send arcs that cross horizon to dummy end time
			a_prime = targettsn.arcs[targettsn.nodes[l1,t1], targettsn.nodes[l2,dummyendtime]]
			push!(newarcset, a_prime)
		else 
			a_prime = targettsn.arcs[targettsn.nodes[l1,t1], targettsn.nodes[l2,t2]]
			push!(newarcset, a_prime)
		end
	end

	return newarcset 

end

#----------------------------------------------------------------------------------------#

function createfragmentsets_online(currstate, hl, ss, sn, lth, drivergroupnum, driversingroup, drivergroupdesc, driverarcs, ghostdriverarcs, odpairs, matchabletrips, enumeratestandardjourneys_flag)
	                                  
	journeys = []

	#Get driver group information
	startloc, starttime = nodesLookup[sn]
	g = drivergroupnum[hl,ss,sn,lth]
	d_ex = driversingroup[hl,ss,sn,lth][1] #Example driver for all the annoying sets that are indexed by driver instead of group
	upcomingfragmentstarttime = [t for t in currstate.T_on_0[d_ex] if (t > starttime) & (t < horizon)]
	
	#----------------------------------------------------------------------#

	starttime_standard = time()
	#Find journeys beginning from driver start node

	if enumeratestandardjourneys_flag == 1
		#Find deadline and round down to beginning of an off period
		if (sn == nodes[hl,0]) & (0 in currstate.T_off[ss])
			deadlinetime = minimum(currstate.T_on_0[d_ex])
		else
			deadlinetime = lth + (shiftlength + 24*maxnightsaway)
		end
		offtimestarttimes = currstate.T_off_0[d_ex] #union(currstate.T_off_0[d_ex], last(currstate.T_off_0[d_ex])+24:24:deadlinetime)
		deadlinetime = [t for t in offtimestarttimes if t <= deadlinetime] == [] ? 0 : maximum([t for t in offtimestarttimes if t <= deadlinetime]) + (24-shiftlength)
		if deadlinetime > horizon
			#Set the deadline at dummy extended horizon (basically, deal with it later)
			deadlinenode = extendednodes[hl, dummyendtime]
		else
			deadlinenode = extendednodes[hl, deadlinetime]
		end

		#Reduce arcs for the driver group based on the deadline
		#timespacenetwork_providetsn(ghosttsn, "outputs/viz/aaaaghost.png", [ghostdriverarcs.A[hl,ss]], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
		availablearcs_ghost = shortestpathonTSN(sn, hl, deadlinetime, ghostdriverarcs.A[hl,ss], ghosttsn)
		#timespacenetwork_providetsn(ghosttsn, "outputs/viz/aaaaghost.png", [availablearcs_ghost], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
		availablearcs_base = translatearcsbetweenTSNs_withextendedarcs(availablearcs_ghost, ghosttsn, basetsn)
		#timespacenetwork("outputs/viz/aaa_all.png", [availablearcs_base], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
		journeyarcs = formatarcset(availablearcs_base, basetsn)

		#Find all journeys beginning from startnode, sn
		journeyendingnodelist = union(N_end, [nodes[hl, t] for t in [t2 for t2 in currstate.T_on_0[d_ex] if t2 <= horizon] if sn < nodes[hl, t] <= deadlinenode], numnodes+1:extendednumnodes)
		if sn < nodes[hl,currstate.T_on_0[d_ex][1]] < journeyendingnodelist[1]
			journeyendingnodelist = union(journeyendingnodelist, nodes[hl,currstate.T_on_0[d_ex][1]])
		end
		journeylist = enumeratealljourneys(sn, journeyendingnodelist, journeyarcs)
		journeys = union(journeys, journeylist)
		
		#include("scripts/onlineimplementation/initializecurrentstatearcs.jl")
		#include("scripts/visualizations/timespacenetwork.jl")
		#timespacenetwork("outputs/viz/aaa_all.png", [availablearcs_base], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
		#for j in 1:20:length(journeylist)
		#	timespacenetwork(string("outputs/viz/aaa_all_",j,".png"), [journeylist[j]], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
		#end

		#----------------------------------------------------------------------#

		#Add journeys starting from future times
		for t in upcomingfragmentstarttime
			fsn = nodes[hl, t] #Future start node
			future_deadlinetime = t + shiftlength + 24*maxnightsaway
			future_deadlinetime = maximum([t for t in offtimestarttimes if t <= future_deadlinetime]) + (24-shiftlength)
			if future_deadlinetime > horizon
				future_deadlinenode = extendednodes[hl, dummyendtime]
			else
				future_deadlinenode = extendednodes[hl, future_deadlinetime]
			end

			#timespacenetwork_providetsn(ghosttsn, "outputs/viz/aaaaghost.png", [ghostdriverarcs.A[hl,ss]], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
			future_availablearcs_ghost = shortestpathonTSN(fsn, hl, future_deadlinetime, ghostdriverarcs.A[hl,ss], ghosttsn)
			#timespacenetwork_providetsn(ghosttsn, "outputs/viz/aaaaghost.png", [future_availablearcs_ghost], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
			future_availablearcs_base = translatearcsbetweenTSNs_withextendedarcs(future_availablearcs_ghost, ghosttsn, basetsn)
			#timespacenetwork("outputs/viz/aaa_all.png", [future_availablearcs_base], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)

			future_journeyarcs = formatarcset(future_availablearcs_base, basetsn)
			future_journeyendingnodelist = union([nodes[hl, t] for t in [t2 for t2 in currstate.T_on_0[d_ex] if t2 <= horizon] if nodes[hl, t] <= future_deadlinenode], N_end, numnodes+1:extendednumnodes)
			future_journeylist = enumeratealljourneys(fsn, future_journeyendingnodelist, future_journeyarcs)
			#for j in 1:20:length(future_journeylist)
			#	timespacenetwork(string("outputs/viz/aaa_all_",j,".png"), [intersect(future_journeylist[j],1:numarcs)], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
			#end
			journeys = union(journeys, future_journeylist)
		end
	else
		
		#Find shortest path home if away
		if (sn == nodes[hl,0]) & (0 in currstate.T_off[ss])
			deadlinetime = minimum(currstate.T_on_0[d_ex])
		else
			deadlinetime = lth + (shiftlength + 24*maxnightsaway)
		end
		offtimestarttimes = currstate.T_off_0[d_ex] #union(currstate.T_off_0[d_ex], last(currstate.T_off_0[d_ex])+24:24:deadlinetime)
		deadlinetime = maximum([t for t in offtimestarttimes if t <= deadlinetime]) + (24-shiftlength)
		if deadlinetime > horizon
			#Set the deadline at dummy extended horizon (basically, deal with it later)
			deadlinenode = extendednodes[hl, dummyendtime]
		else
			deadlinenode = extendednodes[hl, deadlinetime]
		end
		#availablearcs_ghost = shortestpathonTSN(sn, hl, deadlinetime, ghostdriverarcs.A[hl,ss], ghosttsn)
		sphomearcs_ghost, pathendtime = findshortestjourneyhome(sn, hl, deadlinetime, ghostdriverarcs.A[hl,ss], ghosttsn)
		if pathendtime <= horizon
			for t in pathendtime:tstep:minimum(union([t2-tstep for t2 in currstate.T_on_0[d_ex] if pathendtime <= t2 <= horizon], horizon-tstep))
				push!(sphomearcs_ghost, ghosttsn.arcs[ghosttsn.nodes[hl,t],ghosttsn.nodes[hl,t+tstep]])
			end
		else
			throw(DomainError(d_ex, "There is a driver who cannot come home before horizon!"))
		end
		if sphomearcs_ghost != []
			sphomearcs_base = translatearcsbetweenTSNs_withextendedarcs(sphomearcs_ghost, ghosttsn, basetsn)
			#journeyarcs = formatarcset(sphomearcs_base, basetsn)
			push!(journeys, sphomearcs_base)
		end
		
		#Find all idle journeys
		idlejourneys = []

		#Current
		journeyarcs = []
		for t2 in 0:tstep:min(horizon-tstep, 24-tstep)
			push!(journeyarcs, arcs[nodes[hl,t2], nodes[hl,t2+tstep]])	
		end
		push!(idlejourneys, journeyarcs)		

		#Future
		for t in currstate.T_on_0[d_ex]
			journeyarcs = []
			for t2 in t:tstep:min(horizon-tstep, t+24-tstep)
				push!(journeyarcs, arcs[nodes[hl,t2], nodes[hl,t2+tstep]])	
			end
			if journeyarcs != []
				push!(idlejourneys, journeyarcs)
			end
		end

		if idlejourneys != []
			journeys = union(journeys, idlejourneys)
		end
		#println("Before = ", length(journeys))
	end
	time_standard = time() - starttime_standard 
	#----------------------------------------------------------------------#
	
	starttime_extended = time()

	order_journeys, repos_journeys = [], []
	odpairs_homeorig = [(o,d) for (o,d) in odpairs if o==hl]
	odpairs_homedest = [(o,d) for (o,d) in odpairs if d==hl]

	#Driver home location = order origin
	for (o1, d1) in odpairs_homeorig, (o2,d2) in matchabletrips[o1,d1], t1 in [t for t in currstate.T_on_0[d_ex] if t < horizon]
		if (arcLength[o1,d1] > shiftlength) || (arcLength[o2,d2] > shiftlength)

			#Arc locations 
			l1, l2, l3, l4, l5, l6, l7, l8, l9 = o1, d1, d1, o2, o2, d2, d2, o1, o1
			journeylocations = [l1, l2, l3, l4, l5, l6, l7, l8, l9]

			#Calculate arc times
			#Note that this relies on the assumption that the repositioning trips take at most 1 tstep to complete!!!
			t2 = t1 + arcLength[o1, d1] <= horizon ? t1 + arcLength[l1,l2] : dummyendtime #o1 --> d1
			upcomingworkinghours = [t_prime for t_prime in setdiff(0:tstep:horizon, currstate.T_off[ss]) if t_prime >= t2]
			t3 = minimum(union(upcomingworkinghours, dummyendtime)) #Find start time of d1 --> o2
			t4 = d1 == o2 ? t3 : t3 + arcLength[d1,o2] <= horizon ? t3 + arcLength[d1,o2] : dummyendtime #d1 --> o2
			upcomingshifts = [t_prime for t_prime in currstate.T_on_0[d_ex] if t_prime >= t4]
			t5 = minimum(union(upcomingshifts, dummyendtime)) <= horizon ? minimum(union(upcomingshifts, dummyendtime)) : dummyendtime
			t6 = t5 + arcLength[o2,d2] <= horizon ? t5 + arcLength[o2,d2] : dummyendtime #o2 --> d2
			upcomingworkinghours = [t_prime for t_prime in setdiff(0:tstep:horizon, currstate.T_off[ss]) if t_prime >= t6]
			t7 = minimum(union(upcomingworkinghours, dummyendtime)) #Find start time of d2 --> o1
			t8 = d2 == o1 ? t7 : t7 + arcLength[d2,o1] <= horizon ?  t7 + arcLength[d2,o1] : dummyendtime #d2 --> o1
			upcomingshifts = [t_prime for t_prime in currstate.T_on_0[d_ex] if t_prime >= t8]
			t9 = minimum(union(upcomingshifts, dummyendtime)) <= horizon ? minimum(union(upcomingshifts, dummyendtime)) : dummyendtime

			#Create journey
			journeytimes = [t1,t2,t3,t4,t5,t6,t7,t8,t9]
			finaltimeindex = argmax([t for t in journeytimes]) #Stop times at first dummyendtime
			if (dummyendtime in journeytimes) & (finaltimeindex > 1) && (journeylocations[finaltimeindex] == journeylocations[finaltimeindex-1]) 
				#If the last arc before the time horizon ends would be an idle arc, then just drop it
				finaltimeindex = finaltimeindex - 1
			end
			#Nodes
			journeypairlist = [(journeylocations[t], journeytimes[t]) for t in 1:finaltimeindex]
			journeynodelist = [extendednodes[journeypairlist[1][1],journeypairlist[1][2]]]
			for nindex in 2:finaltimeindex
				loc, tm = journeypairlist[nindex]
				lastloc, lasttm = journeypairlist[nindex-1]
				if lastloc == loc
					for t in lasttm:tstep:tm
						push!(journeynodelist, extendednodes[loc,t])
					end
				else
					push!(journeynodelist, extendednodes[loc,tm])
				end
			end
			journeynodelist = unique(journeynodelist)
			#Arcs
			journeyarclist = []
			for n in 1:length(journeynodelist)-1
				push!(journeyarclist, extendedarcs[journeynodelist[n], journeynodelist[n+1]])
			end

			#timespacenetwork(string("outputs/viz/aaa_order.png"), [intersect(journeyarclist,1:numarcs)], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
			push!(repos_journeys, journeyarclist)
		end
	end

	#Driver home location = order destination
	for (o2,d2) in odpairs_homedest, (o1,d1) in matchabletrips[o2,d2], t1 in [t for t in currstate.T_on_0[d_ex] if t < horizon]
		if (arcLength[o1,d1] > shiftlength) || (arcLength[o2,d2] > shiftlength)

			#Arc locations 
			l1, l2, l3, l4, l5, l6, l7, l8, l9 = d2, o1, o1, d1, d1, o2, o2, d2, d2
			journeylocations = [l1, l2, l3, l4, l5, l6, l7, l8, l9]

			#Calculate arc times
			#Note that this relies on the assumption that the repositioning trips take at most 1 tstep to complete!!!
			t2 = d2 == o1 ? t1 : t1 + arcLength[d2, o1] <= horizon ? t1 + arcLength[d2, o1] : dummyendtime #d2 --> o1
			upcomingshifts = [t_prime for t_prime in currstate.T_on_0[d_ex] if t_prime >= t2]
			t3 = minimum(union(upcomingshifts, dummyendtime)) <= horizon ? minimum(union(upcomingshifts, dummyendtime)) : dummyendtime
			t4 = t3 + arcLength[o1,d1] <= horizon ? t3 + arcLength[o1,d1] : dummyendtime #o1 --> d1
			upcomingworkinghours = [t_prime for t_prime in setdiff(0:tstep:horizon, currstate.T_off[ss]) if t_prime >= t4]
			t5 = minimum(union(upcomingworkinghours, dummyendtime)) #Find start time of d1 --> o2
			t6 = d1 == o2 ? t5 : t5 + arcLength[d1,o2] <= horizon ? t5 + arcLength[d1,o2] : dummyendtime #d1 --> o2
			upcomingshifts = [t_prime for t_prime in currstate.T_on_0[d_ex] if t_prime >= t6]
			t7 = minimum(union(upcomingshifts, dummyendtime)) <= horizon ? minimum(union(upcomingshifts, dummyendtime)) : dummyendtime
			t8 = t7 + arcLength[o2,d2] <= horizon ?  t7 + arcLength[o2,d2] : dummyendtime #o2 --> d2
			upcomingshifts = [t_prime for t_prime in currstate.T_on_0[d_ex] if t_prime >= t8]
			t9 = minimum(union(upcomingshifts, dummyendtime)) <= horizon ? minimum(union(upcomingshifts, dummyendtime)) : dummyendtime

			#Create journey
			journeytimes = [t1,t2,t3,t4,t5,t6,t7,t8,t9]
			finaltimeindex = argmax([t for t in journeytimes]) #Stop times at first dummyendtime
			if (dummyendtime in journeytimes) & (finaltimeindex > 1) && (journeylocations[finaltimeindex] == journeylocations[finaltimeindex-1])
				#If the last arc before the time horizon ends would be an idle arc, then just drop it
				finaltimeindex = finaltimeindex - 1
			end
			#Nodes
			journeypairlist = [(journeylocations[t], journeytimes[t]) for t in 1:finaltimeindex]
			journeynodelist = [extendednodes[journeypairlist[1][1],journeypairlist[1][2]]]
			for nindex in 2:finaltimeindex
				loc, tm = journeypairlist[nindex]
				lastloc, lasttm = journeypairlist[nindex-1]
				if lastloc == loc
					for t in lasttm:tstep:tm
						push!(journeynodelist, extendednodes[loc,t])
					end
				else
					push!(journeynodelist, extendednodes[loc,tm])
				end
			end
			journeynodelist = unique(journeynodelist)
			#Arcs
			journeyarclist = []
			for n in 1:length(journeynodelist)-1
				push!(journeyarclist, extendedarcs[journeynodelist[n], journeynodelist[n+1]])
			end

			#timespacenetwork(string("outputs/viz/aaa_order.png"), [intersect(journeyarclist,1:numarcs)], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)

			push!(repos_journeys, journeyarclist)
		end
	end
	journeys = union(journeys, repos_journeys)

	#Stranded driver needs a ride home
	stranded_journeys = []
	if (startloc != hl) # & (arcLength[startloc,hl] > shiftlength)
		o1, d1 = hl, startloc
		backhometrips = [(d1,o1)]
		backhometrips = union(backhometrips, [(o2,d2) for (o2,d2) in odpairs if (distbetweenlocs[o1,d2] <= maxrepositioningdistance) & (distbetweenlocs[o2,d1] <= maxrepositioningdistance) & (arcLength[o1,d2] <= tstep) & (arcLength[o2,d1] <= tstep)])
		#backhometrips = []
		#try 
		#	backhometrips = union(backhometrips, matchabletrips[o1,d1])
		#catch 
		#	backhometrips = union(backhometrips, [(startloc, hl)])
		#and
		for (o2, d2) in backhometrips

			t4 = starttime

			#Arc locations 
			l1, l2, l3, l4, l5, l6, l7, l8, l9 = d2, o1, o1, d1, d1, o2, o2, d2, d2
			journeylocations = [l4, l5, l6, l7, l8, l9]

			#Calculate arc times
			#Note that this relies on the assumption that the repositioning trips take at most 1 tstep to complete!!!
			upcomingworkinghours = [t_prime for t_prime in setdiff(0:tstep:horizon, currstate.T_off[ss]) if t_prime >= t4]
			t5 = minimum(union(upcomingworkinghours, dummyendtime)) #Find start time of d1 --> o2
			t6 = d1 == o2 ? t5 : t5 + arcLength[d1,o2] <= horizon ? t5 + arcLength[d1,o2] : dummyendtime #d1 --> o2
			upcomingshifts = [t_prime for t_prime in currstate.T_on_0[d_ex] if t_prime >= t6]
			t7 = minimum(union(upcomingshifts, dummyendtime)) <= horizon ? minimum(union(upcomingshifts, dummyendtime)) : dummyendtime
			t8 = t7 + arcLength[o2,d2] <= horizon ?  t7 + arcLength[o2,d2] : dummyendtime #o2 --> d2
			upcomingshifts = [t_prime for t_prime in currstate.T_on_0[d_ex] if t_prime >= t8]
			t9 = minimum(union(upcomingshifts, dummyendtime)) <= horizon ? minimum(union(upcomingshifts, dummyendtime)) : dummyendtime

			#Create journey
			journeytimes = [t4,t5,t6,t7,t8,t9]
			finaltimeindex = argmax([t for t in journeytimes]) #Stop times at first dummyendtime
			if (dummyendtime in journeytimes) & (finaltimeindex > 1) && (journeylocations[finaltimeindex] == journeylocations[finaltimeindex-1])
				#If the last arc before the time horizon ends would be an idle arc, then just drop it
				finaltimeindex = finaltimeindex - 1
			end
			#Nodes
			journeypairlist = [(journeylocations[t], journeytimes[t]) for t in 1:finaltimeindex]
			journeynodelist = [extendednodes[journeypairlist[1][1],journeypairlist[1][2]]]
			for nindex in 2:finaltimeindex
				loc, tm = journeypairlist[nindex]
				lastloc, lasttm = journeypairlist[nindex-1]
				if lastloc == loc
					for t in lasttm:tstep:tm
						push!(journeynodelist, extendednodes[loc,t])
					end
				else
					push!(journeynodelist, extendednodes[loc,tm])
				end
			end
			journeynodelist = unique(journeynodelist)
			#Arcs
			journeyarclist = []
			for n in 1:length(journeynodelist)-1
				push!(journeyarclist, extendedarcs[journeynodelist[n], journeynodelist[n+1]])
			end

			#timespacenetwork(string("outputs/viz/aaa_order.png"), [intersect(journeyarclist,1:numarcs)], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)

			push!(stranded_journeys, journeyarclist)
		end
	end
	journeys = union(journeys, stranded_journeys)

	#println("Before = ", length(journeys))
		
	#=order_journeys = []
	odpairs_relevant = union([(o,d) for (o,d) in odpairs if o==hl], [(o,d) for (o,d) in odpairs if d==hl])
	for (o1, d1) in odpairs_relevant, t1 in [t for t in currstate.T_on_0[ss] if t < horizon]

		l1 = hl
		l2 = setdiff([o1, d1], hl)[1]
		#println("$l1, $l2")

		#Calculate arc times
		t2 = t1 + arcLength[l1,l2] <= horizon ? t1 + arcLength[l1,l2] : dummyendtime
		t3 = minimum(union(dummyendtime, [t for t in currstate.T_on_0[ss] if (t >= t2) & (t <= horizon)]))
		t4 = t3 + arcLength[l2,l1] <= horizon ? t3 + arcLength[l2,l1] : dummyendtime
		t5 = minimum(union(dummyendtime, [t for t in currstate.T_on_0[ss] if (t >= t4) & (t <= horizon)]))

		#Find the arcs
		orderarc = extendedarcs[extendednodes[l1,t1], extendednodes[l2,t2]]
		if t2 >= horizon 
			middlearcs, returnarc, endarcs = [], [], []
		elseif t3 >= horizon
			middlearcs = [extendedarcs[extendednodes[l2,t], extendednodes[l2,t+tstep]] for t in t2:tstep:min(t3-tstep, horizon-tstep)]
			returnarc, endarcs = [], []
		elseif t4 >= horizon
			middlearcs = [extendedarcs[extendednodes[l2,t], extendednodes[l2,t+tstep]] for t in t2:tstep:min(t3-tstep, horizon-tstep)]
			returnarc = [extendedarcs[extendednodes[l2,t3], extendednodes[l1,t4]]]
			endarcs = []
		else
			middlearcs = [extendedarcs[extendednodes[l2,t], extendednodes[l2,t+tstep]] for t in t2:tstep:min(t3-tstep, horizon-tstep)]
			returnarc = [extendedarcs[extendednodes[l2,t3], extendednodes[l1,t4]]]
			endarcs = [extendedarcs[extendednodes[l1,t], extendednodes[l1,t+tstep]] for t in t4:tstep:min(t5-tstep, horizon-tstep)]
		end

		#Add the arcs to the journey
		order_journeylist = [orderarc]
		for a in middlearcs
			push!(order_journeylist, a)
		end
		for a in returnarc
			push!(order_journeylist, a)
		end
		for a in endarcs
			push!(order_journeylist, a)
		end

		push!(order_journeys, order_journeylist)
	end
	journeys = union(journeys, order_journeys)
	=#
	
	#=for j in 1:15:length(journeys)
		timespacenetwork(string("outputs/viz/aaa_all_",j,".png"), [intersect(journeys[j],1:numarcs)], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
	end
	myarcs = []
	for item in journeys
		myarcs = union(myarcs, item)
	end
	timespacenetwork("outputs/viz/aaa_all.png", [myarcs], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
	=#

	time_extended = time() - starttime_extended 

	#----------------------------------------------------------------------#

	return journeys, time_standard, time_extended

end 

#-----------------------------------------------------------------------------------------#

function findmatchablepointtopointtrips(odpairs)

	matchabletrips = Dict()
	for (o1,d1) in odpairs
		matchabletrips[o1,d1] = []
	end
	for (o1,d1) in odpairs, (o2,d2) in odpairs 
		if (distbetweenlocs[o1,d2] <= maxrepositioningdistance) & (distbetweenlocs[o2,d1] <= maxrepositioningdistance) & (arcLength[o1,d2] <= tstep) & (arcLength[o2,d1] <= tstep)
			push!(matchabletrips[o1,d1], (o2,d2))
		end
	end
	for (o1,d1) in odpairs
		matchabletrips[o1,d1] = push!(matchabletrips[o1,d1], (d1,o1))
	end
	for (o1,d1) in odpairs
		matchabletrips[o1,d1] = unique(matchabletrips[o1,d1])
	end

	return matchabletrips

end

#-----------------------------------------------------------------------------------------#

function initializedriversetjourneys(currstate, driversets, drivergroupnum, driversingroup, drivergroupdesc, driverarcs, ghostdriverarcs, enumeratestandardjourneys_flag)
		
	#(hl,ss,sn,lth) = (18, 1, 5, -36) driversets[1]
	numfragments = Dict()
	fragmentscontaining = Dict()
	fragmentarcs = Dict()
	F_plus_g, F_minus_g = Dict(), Dict()
	N_flow_g = Dict()

    odpairs = unique([(originloc[i], destloc[i]) for i in currstate.orders])
	matchabletrips = findmatchablepointtopointtrips(odpairs)

	totalstandardtime = 0
	totalextendedtime = 0
		
	for (hl,ss,sn,lth) in driversets
		#println("$hl,$ss,$sn,$lth")
        numfragments[hl,ss,sn,lth] = 0
        for a in 1:extendednumarcs
            fragmentscontaining[hl,ss,sn,lth,a] = []
        end
        for n in 1:extendednumnodes
            F_plus_g[hl,ss,sn,lth,n] = []
            F_minus_g[hl,ss,sn,lth,n] = []
        end
    
        N_flow_g[hl,ss,sn,lth] = []
        for n in setdiff(1:extendednumnodes, union(sn, N_end))
            push!(N_flow_g[hl,ss,sn,lth], n)
        end

		journeys, time_standard, time_extended = createfragmentsets_online(currstate, hl,ss,sn,lth, drivergroupnum, driversingroup, drivergroupdesc, driverarcs, ghostdriverarcs, odpairs, matchabletrips, enumeratestandardjourneys_flag)
		totalstandardtime += time_standard
		totalextendedtime += time_extended

        #Process the list of journeys
        for journey in journeys
            
            #Get journey number
            numfragments[hl,ss,sn,lth] += 1
            j = numfragments[hl,ss,sn,lth]
            
            #Find arcs contained on journey
			fragmentarcs[hl,ss,sn,lth,j] = journey
			for a in journey
				push!(fragmentscontaining[hl,ss,sn,lth,a], j)
			end

            #Get journey flow balance
			n1 = arcLookup[first(journey)][1]
			n2 = arcLookup[last(journey)][2]
            push!(F_plus_g[hl,ss,sn,lth,n1], j)
            push!(F_minus_g[hl,ss,sn,lth,n2], j)
        end

	end

	println("Time to enumerate standard journeys = ", totalstandardtime)
	println("Time to enumerate extended journeys = ", totalextendedtime)

    return numfragments, fragmentscontaining, fragmentarcs, F_plus_g, F_minus_g, N_flow_g

end

#-----------------------------------------------------------------------------------------#

function getfragmentstats(currstate, driversets, numfragments, fragmentarcs)

    fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = Dict(), Dict(), Dict(), Dict()

	for (hl,ss,sn,lth) in driversets
        workingfragments[hl,ss,sn,lth] = []
	end
	
	for (hl,ss,sn,lth) in driversets, f in 1:numfragments[hl,ss,sn,lth] 
        
        drivinghours, workinghours = 0, 0
		fragmentstarttime = nodesLookup[arcLookup[fragmentarcs[hl,ss,sn,lth,f][1]][1]][2]
		fragmentendtime = min(horizon, nodesLookup[arcLookup[last(fragmentarcs[hl,ss,sn,lth,f])][2]][2])
        for a in fragmentarcs[hl,ss,sn,lth,f]
            l1,t1 = nodesLookup[arcLookup[a][1]]
            l2,t2 = nodesLookup[arcLookup[a][2]]
            if !(t1 in currstate.T_off[ss]) & (l1 != l2)
                drivinghours += t2-t1
                workinghours += t2-t1
            elseif !(t1 in currstate.T_off[ss]) & (l1 == l2) & !(l1 == hl)
                workinghours += t2-t1
            end
			#Check me!!
			fragmentnightsaway[hl,ss,sn,lth,f] = ceil((fragmentendtime - fragmentstarttime) / 24) - 1
			#if (t1 in currstate.T_off_0[ss]) & (l1 != hl)
			#	fragmentnightsaway[hl,ss,sn,lth,f] += 1
			#end
        end

        fragdrivinghours[hl,ss,sn,lth,f] = drivinghours
        fragworkinghours[hl,ss,sn,lth,f] = workinghours

        if workinghours > 1e-4
            push!(workingfragments[hl,ss,sn,lth], f)
        end
	end

    return fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway

end

#-----------------------------------------------------------------------------------------#

function initializecurrentstatearcs(currstate, enumeratestandardjourneys_flag)

	#include("scripts/onlineimplementation/initializecurrentstatearcs.jl")

	println("Primary")
	@time primaryarcs, extendedtimearcs, orderarcs, driverarcs, hasdriverarcs, ghostdriverarcs = initializearcsets(A_space, A_plus, A_minus, currstate.orders, currstate.Origin, currstate.Destination, currstate.driverStartNodes, currstate.T_off)
	#R_off = findreturnhomearcsets(driverarcs, currstate.T_off_constr)
	println("MAG")
	if operations == "relay"
	    @time magarcs = initializeorderarcsets(k, currstate.orders, originloc, destloc, currstate.Origin, currstate.Destination, currstate.shortesttriptimes)
    elseif operations == "ptp"
		@time magarcs = orderarcs
	end
	#driversets, driverSetStartNodes, numfragments, fragmentscontaining, F_plus_ls, F_minus_ls, N_flow_ls, numeffshifts, effshift, shiftsincluded, fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = initializejourneymodel(maxnightsaway, currstate.T_off, currstate.T_on_0)
    println("Driver sets")
    #Create driver sets and journeys
    @time driversets, driversingroup, numdrivergroups, drivergroupnum, drivergroupdesc, numeffshifts, effshift, shiftsincluded = finddriversets_online(currstate.T_off, currstate.driverStartNodes, currstate.lasttimehome) 
	journeystart = time()
	numfragments, fragmentscontaining, fragmentarcs, F_plus_g, F_minus_g, N_flow_g = initializedriversetjourneys(currstate, driversets, drivergroupnum, driversingroup, drivergroupdesc, driverarcs, ghostdriverarcs, enumeratestandardjourneys_flag)
	journeytime = time() - journeystart

	#=

maxnightsaway = 3
maxrepositioningdistance = 100

driversets, driversingroup, numdrivergroups, drivergroupnum, drivergroupdesc, numeffshifts, effshift, shiftsincluded = finddriversets_online(currstate.T_off, currstate.driverStartNodes, currstate.lasttimehome) 

include("scripts/onlineimplementation/initializecurrentstatearcs.jl")

journeystart = time()
numfragments, fragmentscontaining, fragmentarcs, F_plus_g, F_minus_g, N_flow_g = initializedriversetjourneys(currstate, driversets, drivergroupnum, driversingroup, drivergroupdesc, driverarcs, ghostdriverarcs, enumeratestandardjourneys_flag)
journeytime = time() - journeystart

fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = getfragmentstats(currstate, driversets, numfragments, fragmentarcs)
currarcs = (orderarcs=orderarcs, driverarcs=driverarcs, hasdriverarcs=hasdriverarcs, magarcs=magarcs); #, R_off=R_off)
currfragments = (driversets=driversets, driversingroup=driversingroup, numdrivergroups=numdrivergroups, drivergroupnum=drivergroupnum, drivergroupdesc=drivergroupdesc, numeffshifts=numeffshifts, effshift=effshift, shiftsincluded=shiftsincluded,numfragments=numfragments, fragmentscontaining=fragmentscontaining, fragmentarcs=fragmentarcs, F_plus_g=F_plus_g, F_minus_g=F_minus_g, N_flow_g=N_flow_g, fragdrivinghours=fragdrivinghours, fragworkinghours=fragworkinghours, workingfragments=workingfragments, fragmentnightsaway=fragmentnightsaway);
totalnumjourneys = sum(currfragments.numfragments[ds] for ds in currfragments.driversets)


include("scripts/journeybasedmodel/solvejourneymodel_online.jl")		
ip_obj, x_ip, z_ip, ip_time, ip_bound = solvejourneymodel(0, opt_gap, currstate, currarcs, currfragments, orderarcs)
#ip_time = 0

totalnightsaway = 0
for (i1,i2,i3,i4) in currfragments.driversets, f in 1:currfragments.numfragments[i1,i2,i3,i4]
	if value(z_ip[(i1,i2,i3,i4), f]) > 1e-4
		totalnightsaway += fragmentnightsaway[i1,i2,i3,i4,f]
	end
end
totalnightsaway / length(drivers)

println("$ex, $maxnightsaway, $totalnumjourneys, $journeytime")
df = DataFrame(ex=ex, maxnightsaway=maxnightsaway, journeys=totalnumjourneys, calctime=journeytime, solvetime=ip_time)
CSV.write("outputs/viz/ptpjourneyexamples/journeystats.csv", df, append=true)

=#

	@time fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = getfragmentstats(currstate, driversets, numfragments, fragmentarcs)

	fragments = Dict()
	for item in driversets
		fragments[item] = [j for j in 1:numfragments[item]]
	end
		
    currarcs = (orderarcs=orderarcs, driverarcs=driverarcs, hasdriverarcs=hasdriverarcs, magarcs=magarcs, ghostdriverarcs=ghostdriverarcs) #, R_off=R_off)
    #currfragments = (driversets=driversets, driverSetStartNodes=driverSetStartNodes, numfragments=numfragments, 
    #            fragmentscontaining=fragmentscontaining, F_plus_ls=F_plus_ls, F_minus_ls=F_minus_ls, N_flow_ls=N_flow_ls, 
    #            effshift=effshift, shiftsincluded=shiftsincluded, fragdrivinghours=fragdrivinghours, 
    #            fragworkinghours=fragworkinghours, workingfragments=workingfragments, fragmentnightsaway=fragmentnightsaway)
    currfragments = (driversets=driversets, driversingroup=driversingroup, numdrivergroups=numdrivergroups, drivergroupnum=drivergroupnum, 
                drivergroupdesc=drivergroupdesc, numeffshifts=numeffshifts, effshift=effshift, shiftsincluded=shiftsincluded,
                fragments=fragments, numfragments=numfragments, fragmentscontaining=fragmentscontaining, fragmentarcs=fragmentarcs, 
                F_plus_g=F_plus_g, F_minus_g=F_minus_g, N_flow_g=N_flow_g, 
                #effshift=effshift, shiftsincluded=shiftsincluded, 
                fragdrivinghours=fragdrivinghours, fragworkinghours=fragworkinghours, 
                workingfragments=workingfragments, fragmentnightsaway=fragmentnightsaway)

    return currarcs, currfragments, primaryarcs, extendedtimearcs, numeffshifts, journeytime

end

#=
2,2,2,0
(hl,ss,sn,lth) = driversets[14]
myarcs = []
for j in 1:currfragments.numfragments[(hl,ss,sn,lth)]
	myarcs = union(myarcs, currfragments.fragmentarcs[hl,ss,sn,lth,j])
end
include("scripts/visualizations/timespacenetwork.jl")
timespacenetwork("outputs/viz/aaa_all.png", [myarcs], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
for j in 1:10000:currfragments.numfragments[10, 2, 10, 0]
	timespacenetwork(string("outputs/viz/aaa_all_",j,".png"), [ currfragments.fragmentarcs[10, 2, 10, 0,j]], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
end
=#

