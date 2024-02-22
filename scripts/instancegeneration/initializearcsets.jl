
function orderarcreduction(orders, Origin, Destination)
	
	#For now, we are using shortest path regardless of driver availability (this would be more time consuming to do)
	traveltime_rdd = cacheShortestTravelTimes(numlocs, prearcs, "rdd time")
	traveltime_llr = cacheShortestTravelTimes(numlocs, prearcs, "llr time")
	orderarcset, orderarcset_space, A_plus_i, A_minus_i = Dict(), Dict(), Dict(), Dict()

	for i in orders
		orderarcset[i] = [dummyarc]
		orderarcset_space[i] = []
	end
	for i in orders, n in 1:extendednumnodes
		if n == Origin[i][1]
			A_plus_i[i,n] = [dummyarc]
			A_minus_i[i,n] = []
		elseif n == last(Destination[i])
			A_plus_i[i,n] = []
			A_minus_i[i,n] = [dummyarc]
		else 
			A_plus_i[i,n] = []
			A_minus_i[i,n] = []
		end
	end

	for i in orders
		destinationlocation = nodesLookup[Destination[i][1]][1]
		#push!(Destination[i], extendednodes[destinationlocation, dummyendtime])
		A_minus_i[i, extendednodes[destinationlocation, dummyendtime]] = []
		for n2 in setdiff(N_end, Destination[i])
			arc_ext = extendedarcs[n2, extendednodes[destinationlocation, dummyendtime]]
			push!(orderarcset[i], arc_ext)
			push!(A_plus_i[i, n2], arc_ext)
			push!(A_minus_i[i, extendednodes[destinationlocation, dummyendtime]], arc_ext)
		end
	end

	prearcs_aug = deepcopy(prearcs)
	for l in 1:numlocs
		push!(prearcs_aug, (l, l, tstep, tstep))
	end

	for i in orders, arc in prearcs_aug
		arcorig, arcdest, arctraveltime = arc[1], arc[2], arc[3]
		if (arcdest != orderOriginalStartLoc[i]) & (arcorig != nodesLookup[Destination[i][1]][1]) & (arctraveltime <= shiftlength)
			startloc, starttime = nodesLookup[Origin[i][1]]
			endloc = nodesLookup[Destination[i][1]][1]

			t1 = traveltime_rdd[startloc, arcorig]
			t2 = arc[3]

			for t in 0:tstep:horizon-t2
				if starttime + t1 <= t
					push!(orderarcset[i], arcs[nodes[arcorig, t], nodes[arcdest, t + t2]])
					if arcorig != arcdest
						push!(orderarcset_space[i], arcs[nodes[arcorig, t], nodes[arcdest, t + t2]])
					end				
				end
			end
		end
	end

	#Create A_plus and A_minus lists
	for i in orders, n in 1:numnodes, a in A_plus[n]
		if (a in orderarcset[i]) & !(a in A_plus_i[i,n])
			push!(A_plus_i[i,n], a)
		end
	end
	for i in orders, n in 1:numnodes, a in A_minus[n]
		if (a in orderarcset[i]) & !(a in A_minus_i[i,n])
			push!(A_minus_i[i,n], a)
		end
	end

	for i in orders
		orderarcset[i] = union(sort(setdiff(orderarcset[i], dummyarc), by = x -> nodesLookup[arcLookup[x][1]][2]), dummyarc)
	end

	return orderarcset, orderarcset_space, A_plus_i, A_minus_i

end

#---------------------------------------------------------------------------------------#

function driverarcreduction(driverStartNodes, T_off)

	traveltime = cacheShortestTravelTimes(numlocs, prearcs, "rdd time")
	homeArcSet, homeArcSet_space, availableDrivers, A_plus_d, A_minus_d = Dict(), Dict(), Dict(), Dict(), Dict()
	fullhourlist = [t for t in 0:tstep:horizon-tstep]

	for d in drivers
		homeArcSet[d] = []
		homeArcSet_space[d] = []
	end
	for a in 1:numarcs
		availableDrivers[a] = []
	end
	for d in drivers, n in 1:numnodes
		A_plus_d[d,n] = []
		A_minus_d[d,n] = []
	end

	prearcs_aug = deepcopy(prearcs)
	for l in 1:numlocs
		push!(prearcs_aug, (l, l, tstep, tstep))
	end

	closelocs = Dict()
	for d in drivers
		closelocs[d] = []
	end

	for d in drivers, arc in prearcs_aug
		#Travel times between all relevant locations
		orig, dest = arc[1], arc[2]
		h = driverHomeLocs[d]
		t1 = traveltime[h, orig]
		t2 = arc[3]
		t3 = traveltime[dest, h]

		#Check earliest we can get to the arc origin location from our current location
		driverstartingloc, driverstartingtime = nodesLookup[driverStartNodes[d]]
		firstshifttime = setdiff(fullhourlist, T_off[drivershift[d]])[1]
		drivinghourstoarcorig = traveltime[driverstartingloc, orig]
		totalhourstoarcorig = floor(drivinghourstoarcorig/shiftlength) * (24 - shiftlength) + drivinghourstoarcorig
		earliesttime = totalhourstoarcorig + max(firstshifttime, driverstartingtime)

		#If the arc is feasible for the driver, add it to the driver's arc list
		if ((t1 + t2 <= shiftlength) & (t3 <= shiftlength)) || ((t1 <= shiftlength) & (t2 + t3 <= shiftlength))
			if !(orig in closelocs[d])
				push!(closelocs[d], orig)
			end
			if !(dest in closelocs[d])
				push!(closelocs[d], dest)
			end

			for t in setdiff([t4 for t4 in fullhourlist if t4 >= earliesttime], T_off[drivershift[d]])
				#Arc must finish before the end of the horizon and before the "next" off hour of the driver
				if t + t2 <= [t4 for t4 in union(T_off[drivershift[d]], horizon) if t4 > t][1]
					startnode, endnode = nodes[orig, t],  nodes[dest, t + t2]
					a = arcs[startnode, endnode]
					push!(homeArcSet[d], a)
					if orig != dest
						push!(homeArcSet_space[d], a)
					end				
					push!(availableDrivers[a], d)
					push!(A_plus_d[d, startnode], a)
					push!(A_minus_d[d,endnode], a)
				end
			end
		end
	end
	
	#Add stay-at-home arcs for driver's off hours
	for d in drivers, t in setdiff(T_off[drivershift[d]], [horizon]), l in closelocs[d]
		startnode, endnode = nodes[l, t], nodes[l, t + tstep]
		a = arcs[startnode, endnode]
		push!(homeArcSet[d], a)
		push!(availableDrivers[a], d)
		push!(A_plus_d[d, startnode], a)
		push!(A_minus_d[d,endnode], a)
	end

	return homeArcSet, homeArcSet_space, availableDrivers, A_plus_d, A_minus_d, closelocs

end

#---------------------------------------------------------------------------------------#

function truckarcreduction(A_space, A_plus, A_minus, availabledrivers)

	A_hasdriver = []
	yupperbound = []
	A_hasdriver_space = []
	for a in 1:numarcs
		ub = 99999
		if a in A_space
			ub = length(availabledrivers[a])
		end
		push!(yupperbound, ub)
		if (a in A_space) & (ub > 0)
			push!(A_hasdriver, a)
			push!(A_hasdriver_space, a)
		elseif !(a in A_space)
			push!(A_hasdriver, a)
		end
	end

	A_plus_hd, A_minus_hd = Dict(), Dict()
	for n in 1:numnodes
		A_plus_hd[n] = []
		A_minus_hd[n] = []
		for a in A_plus[n]
			if a in A_hasdriver
				push!(A_plus_hd[n], a)
			end
		end
		for a in A_minus[n]
			if a in A_hasdriver
				push!(A_minus_hd[n], a)
			end
		end
	end

	return A_hasdriver, yupperbound, A_hasdriver_space, A_plus_hd, A_minus_hd

end

#---------------------------------------------------------------------------------------#

function initializearcsets(A_space, A_plus, A_minus, orders, Origin, Destination, driverStartNodes, T_off)
    
    #Order arc sets - viable paths from origin to destination
    orderarcset_full, orderarcset_space_full, A_plus_i_full, A_minus_i_full = orderarcreduction(orders, Origin, Destination)

    #Identify feasible arcs for each driver and trucks
    if solutionmethod == "nr"
        driverarcset, driverarcset_space, availabledrivers, A_plus_d, A_minus_d, closelocs = driverArcSetsByDriver_nonrelay(numlocs, numarcs, numnodes, prearcs, drivers, tstep, horizon, nodes, arcs, assignedDrivers, A_minus, A_plus, T_off, drivershift, driverHomeLocs, T_off_0, shiftlength)
        A_hasdriver, yupperbound, A_hasdriver_space, A_plus_hd, A_minus_hd = truckarcreduction(A_space, A_plus, A_minus, availabledrivers)
    else
        driverarcset, driverarcset_space, availabledrivers, A_plus_d, A_minus_d, closelocs = driverarcreduction(driverStartNodes, T_off)
        A_hasdriver, yupperbound, A_hasdriver_space, A_plus_hd, A_minus_hd = truckarcreduction(A_space, A_plus, A_minus, availabledrivers)
    end

    #Format arc sets
    primaryarcs = (A=1:numarcs, A_space=A_space, A_minus=A_minus, A_plus=A_plus, available=[], closelocs=[]);
    extendedtimearcs = (A=1:extendednumarcs, A_space=A_space, A_minus=A_minus, A_plus=A_plus, available=[], closelocs=[]);
    orderarcs = (A=orderarcset_full, A_space=orderarcset_space_full, A_minus=A_minus_i_full, A_plus=A_plus_i_full, available=[], closelocs=[]);
    driverarcs = (A=driverarcset, A_space=driverarcset_space, A_minus=A_minus_d, A_plus=A_plus_d, available=availabledrivers, closelocs=closelocs);
    hasdriverarcs = (A=A_hasdriver, A_space=A_hasdriver_space, A_minus=A_minus_hd, A_plus=A_plus_hd, available=[], closelocs=[]);

    return primaryarcs, extendedtimearcs, orderarcs, driverarcs, hasdriverarcs

end




