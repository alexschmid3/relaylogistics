
function findtraveltimesanddistances(orders, Origin, Destination)
	
	#Calculate shortest paths in miles between all pairs of locations
	distbetweenlocs, shortestpatharclists = cacheShortestDistance(numlocs, prearcs)
	traveltimebetweenlocs_rdd = cacheShortestTravelTimes(numlocs, prearcs, "rdd time")
	traveltimebetweenlocs_raw = cacheShortestTravelTimes(numlocs, prearcs, "raw time")
	traveltimebetweenlocs_llr = cacheShortestTravelTimes(numlocs, prearcs, "llr time")

	#Calculate the shortest trip times for each order i 
	shortesttriptimes = []
	for i in orders
		#Using shortest path distances, ignoring driver availability
		if traveltimefordelay_flag == 0
			shortestpathtime = traveltimebetweenlocs_rdd[nodesLookup[Origin[i][1]][1], nodesLookup[Destination[i][1]][1]]
		elseif traveltimefordelay_flag == 1
			shortestpathtime = traveltimebetweenlocs_raw[nodesLookup[Origin[i][1]][1], nodesLookup[Destination[i][1]][1]]
		elseif traveltimefordelay_flag == 2
			shortestpathtime = traveltimebetweenlocs_llr[nodesLookup[Origin[i][1]][1], nodesLookup[Destination[i][1]][1]]
		end
		
		push!(shortesttriptimes, shortestpathtime)
	end

	shortesttripdists = []
	for i in orders
		shortestpathdist = distbetweenlocs[nodesLookup[Origin[i][1]][1], nodesLookup[Destination[i][1]][1]]
		push!(shortesttripdists, shortestpathdist)
	end

	return distbetweenlocs, shortesttriptimes, shortestpatharclists, traveltimebetweenlocs_rdd, traveltimebetweenlocs_raw, traveltimebetweenlocs_llr

end

#---------------------------------------------------------------------------------------#

function extendtimespacenetwork(nodesLookup, arcLookup, A_minus, A_plus, c, u, distbetweenlocs)

	#Add final legs for orders that are unfinished at the end of the horizon
	extendednodes, extendednumnodes, extendedarcs, extendednumarcs = copy(nodes), copy(numnodes), copy(arcs), copy(numarcs)
	for l in 1:numlocs
		nodesLookup[extendednumnodes + 1] = (l, dummyendtime)
		extendednodes[l, dummyendtime] = extendednumnodes + 1
		A_minus[extendednumnodes + 1], A_plus[extendednumnodes + 1] = [], []
		extendednumnodes += 1
	end
	if operations == "ptp"
		for l1 in 1:numlocs, l2 in 1:numlocs, t in max(0,horizon-arcLength[l1,l2]+tstep):tstep:horizon
			n1, n2 = extendednodes[l1, t], extendednodes[l2, dummyendtime]
			extendedarcs[n1,n2] = extendednumarcs + 1
			arcLookup[extendednumarcs + 1] = (n1, n2)
			push!(c, distbetweenlocs[l1,l2] * (1 + finallegdistancepenalty))
			push!(u, taxicostpct * distbetweenlocs[l1,l2] * (1 + finallegdistancepenalty))
			push!(A_minus[n2], extendednumarcs + 1)
			push!(A_plus[n1], extendednumarcs + 1)
			extendednumarcs += 1
		end
	elseif operations == "relay"
		for l1 in 1:numlocs, l2 in 1:numlocs
			n1, n2 = extendednodes[l1, horizon], extendednodes[l2, dummyendtime]
			extendedarcs[n1,n2] = extendednumarcs + 1
			arcLookup[extendednumarcs + 1] = (n1, n2)
			push!(c, distbetweenlocs[l1,l2] * (1 + finallegdistancepenalty))
			push!(u, taxicostpct * distbetweenlocs[l1,l2] * (1 + finallegdistancepenalty))
			push!(A_minus[n2], extendednumarcs + 1)
			push!(A_plus[n1], extendednumarcs + 1)
			extendednumarcs += 1
		end
	end

	return nodesLookup, arcLookup, A_minus, A_plus, c, extendednodes, extendednumnodes, extendedarcs, extendednumarcs, u

end

#---------------------------------------------------------------------------------------#

function extendDestination(orders, Destination, extendednodes)

	for i in orders
		destinationlocation = nodesLookup[Destination[i][1]][1]
		push!(Destination[i], extendednodes[destinationlocation, dummyendtime])
	end

	return Destination

end

#---------------------------------------------------------------------------------------#

function findreturnhomearcsets(driverarcs, T_off_constr)

	#Find return home arc sets
	R_off = Dict()
	for d in drivers
		R_off[d] = []
	end
	for d in drivers, t in [t2 for t2 in T_off_constr[d] if t2 <= horizon-24+tstep] 
		h = driverHomeLocs[d]
		if t + 24 == horizon 
			a1 = arcs[nodes[h,t], nodes[h,t+tstep]]
			arcset = union(a1, [a for a in driverarcs.A_minus[d,nodes[h,horizon]]])
		else
			a1, a2 = arcs[nodes[h,t], nodes[h,t+tstep]], arcs[nodes[h,t+24], nodes[h,t+24+tstep]]
			arcset = [a1, a2]
		end
		push!(R_off[d], arcset)
	end

	return R_off

end

#---------------------------------------------------------------------------------------#

function calcarcfinishtimes()

	#Create dummy arc for column generation
	dummyarc = extendednumarcs + 1
	push!(c, 100000000)
	allarcs = extendednumarcs + 1

	#Finish times of arcs
	arcfinishtime = []
	if traveltimefordelay_flag == 0
		for a in 1:numarcs
			push!(arcfinishtime, nodesLookup[arcLookup[a][2]][2])
		end
		for a in numarcs+1:extendednumarcs
			startloc, endloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
			push!(arcfinishtime, horizon + (1 + finallegtimepenalty) * traveltimebetweenlocs_rdd[startloc, endloc] )
		end
	elseif traveltimefordelay_flag >=1
		for a in 1:numarcs
			push!(arcfinishtime, nodesLookup[arcLookup[a][1]][2] + truetraveltime[a])
		end
		for a in numarcs+1:extendednumarcs
			startloc, endloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
			push!(arcfinishtime, horizon + (1 + finallegtimepenalty) * traveltimebetweenlocs_llr[startloc, endloc] )
		end
	end
	push!(arcfinishtime, 1000)

	arcLookup[dummyarc] = (-1,-2)
	nodesLookup[-1] = (-1,0)
	nodesLookup[-2] = (-1,100000)

	return arcLookup, nodesLookup, arcfinishtime, dummyarc, allarcs

end 

#---------------------------------------------------------------------------------------#

function getdriverandshiftinfo()

	#Driver information
	driversintransit = []
	drivers, driverStartNodes, driverEndNodes, driverHomeLocs, assignedDrivers = readdrivers("data/pilots.csv", maxdrivers, numlocs, nodes, horizon)
	N_flow_t, N_flow_d = flowbalancenodesets_td(drivers, numnodes, driverStartNodes, N_0, N_end)

	#Driver shifts
	alltimeswithinview = max(2*horizon, horizon + 2*24*maxnightsaway, onlinetimehorizon)
	T_off_Monday8am, T_off, drivershift, T_off_0, T_off_constr, numshifts, T_on_0 = createdrivershifts(drivers, shiftlength, tstep, drivershifttstep, alltimeswithinview)

	return driversintransit, drivers, driverStartNodes, driverEndNodes, driverHomeLocs, assignedDrivers, N_flow_t, N_flow_d, alltimeswithinview, T_off_Monday8am, T_off, drivershift, T_off_0, T_off_constr, numshifts, T_on_0

end

#---------------------------------------------------------------------------------------#

function findorderstartsandtransits(orders, Origin)

	orderOriginalStartTime = Dict()
	for i in orders
		orderOriginalStartTime[i] = nodesLookup[Origin[i][1]][2]
	end

	orderintransit_flag = Dict()
	for i in orders
		orderintransit_flag[i] = 0
	end

	return orderOriginalStartTime, orderintransit_flag

end

#---------------------------------------------------------------------------------------#

#Reserve trucks for the in transit orders
function findtrucksintransit(ordersinprogress, originloc, available)
	placeholder, loctruckcounter, trucksintransit = Dict(), zeros(numlocs), []
	for i in ordersinprogress
		currentloc, availtime = originloc[i][1], available[i][1]
		try
			placeholder[currentloc, availtime] += 1
		catch
			placeholder[currentloc, availtime] = 1
		end
	end
	for item in placeholder
		currentloc, availtime, ttltrucks = item[1][1], item[1][2], item[2]
		push!(trucksintransit, (currentloc, availtime, ttltrucks))
		loctruckcounter[currentloc] += ttltrucks
	end

	return loctruckcounter, trucksintransit
end

#---------------------------------------------------------------------------------------#

#Modify m_0 if needed to accomodate the trucks that are already in transit
function adjust_m_0(m_0, loctruckcounter)

	#Note: THIS IS CHANGING THE INITIAL STATE OF THE INSTANCE
	truckexcess = m_0 - loctruckcounter
	for l in 1:numlocs
		if loctruckcounter[l] > m_0[l]
			addltrucks = loctruckcounter[l] - m_0[l]
			for j in 1:addltrucks
				extraindex = argmax(truckexcess)
				truckexcess[extraindex] -= 1
				m_0[extraindex] -= 1
				m_0[l] += 1
			end
		end
	end

	return m_0

end
