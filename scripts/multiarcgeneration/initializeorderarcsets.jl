
function orderarcs_dummyonly(orders, originloc, destloc, Origin, Destination)
	
	orderArcSet, orderArcSet_space, A_plus_i, A_minus_i = Dict(), Dict(), Dict(), Dict()

	#Create dummy arcs for arc lists
	for i in orders
		onearc = extendedarcs[nodes[originloc[i], horizon], extendednodes[destloc[i], dummyendtime]]
		orderArcSet[i] = [dummyarc, onearc]
		orderArcSet_space[i] = []

		for n in 1:extendednumnodes
			if n == Origin[i][1]
				A_plus_i[i,n] = [dummyarc]
				A_minus_i[i,n] = []
			elseif n == nodes[originloc[i], horizon]
				A_plus_i[i,n] = [onearc]
				A_minus_i[i,n] = []
			elseif n == last(Destination[i])
				A_plus_i[i,n] = []
				A_minus_i[i,n] = [dummyarc, onearc]
			else 
				A_plus_i[i,n] = []
				A_minus_i[i,n] = []
			end
		end
	end
	
	return orderArcSet, orderArcSet_space, A_minus_i, A_plus_i

end

#---------------------------------------------------------------------------------------#

function orderarcs_warmstart(k, orders, originloc, destloc, Origin, Destination, shortesttriptimes)
	
	#For now, we are using shortest path regardless of driver availability (this would be more time consuming to do)
	traveltime_rdd = cacheShortestTravelTimes(numlocs, prearcs, "rdd time")
	traveltime_llr = cacheShortestTravelTimes(numlocs, prearcs, "llr time")
	orderArcSet, orderArcSet_space, A_plus_i, A_minus_i = Dict(), Dict(), Dict(), Dict()

	for i in orders
		onearc = extendedarcs[nodes[originloc[i], horizon], extendednodes[destloc[i], dummyendtime]]
		orderArcSet[i] = [dummyarc, onearc]
		orderArcSet_space[i] = []

		for n in 1:extendednumnodes
			if n == Origin[i][1]
				A_plus_i[i,n] = [dummyarc]
				A_minus_i[i,n] = []
			elseif n == nodes[originloc[i], horizon]
				A_plus_i[i,n] = [onearc]
				A_minus_i[i,n] = []
			elseif n == last(Destination[i])
				A_plus_i[i,n] = []
				A_minus_i[i,n] = [dummyarc, onearc]
			else 
				A_plus_i[i,n] = []
				A_minus_i[i,n] = []
			end
		end
	end

	#for i in orders
	#	destinationlocation = nodesLookup[Destination[i][1]][1]
	#	#push!(Destination[i], extendednodes[destinationlocation, dummyendtime])
	#	A_minus_i[i, extendednodes[destinationlocation, dummyendtime]] = []
	#	for n2 in N_end
	#		arc_ext = extendedarcs[n2, extendednodes[destinationlocation, dummyendtime]]
	#		push!(orderArcSet[i], arc_ext)
	#		push!(A_plus_i[i, n2], arc_ext)
	#		push!(A_minus_i[i, extendednodes[destinationlocation, dummyendtime]], arc_ext)
	#	end
	#end

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
			if traveltimefordelay_flag == 0
				t3 = traveltime_rdd[arcdest, endloc]
			elseif traveltimefordelay_flag >= 1
				t3 = traveltime_llr[arcdest, endloc]
			end

			#May want to toggle between these
			if ktype_flag == "abs"
				deadline = starttime + shortesttriptimes[i] + k
			elseif ktype_flag == "pct"
				deadline = starttime + shortesttriptimes[i] + k * shortesttriptimes[i]
			elseif ktype_flag == "min24"
				deadline = max(starttime + shortesttriptimes[i] + k * shortesttriptimes[i], starttime + 24)
			end

			for t in 0:tstep:horizon-t2
				if (starttime + t1 <= t) & (t + t2 + t3 <= deadline)
					push!(orderArcSet[i], arcs[nodes[arcorig, t], nodes[arcdest, t + t2]])
					if arcorig != arcdest
						push!(orderArcSet_space[i], arcs[nodes[arcorig, t], nodes[arcdest, t + t2]])
					end				
				end
			end
		end

	end

	#Create A_plus and A_minus lists
	for i in orders, n in 1:numnodes, a in A_plus[n]
		if (a in orderArcSet[i]) & !(a in A_plus_i[i,n])
			push!(A_plus_i[i,n], a)
		end
	end
	for i in orders, n in 1:numnodes, a in A_minus[n]
		if (a in orderArcSet[i]) & !(a in A_minus_i[i,n])
			push!(A_minus_i[i,n], a)
		end
	end

	return orderArcSet, orderArcSet_space, A_minus_i, A_plus_i

end

#---------------------------------------------------------------------------------------#

#Very similar to the function above so I will use that one for now
#orderarcreduction_unik?
function noideawhatthisishonestly(k, ktype_flag, prearcs, shortesttriptimes)
	
	#For now, we are using shortest path regardless of driver availability (this would be more time consuming to do)
	traveltime_rdd = cacheShortestTravelTimes(numlocs, prearcs, "rdd time")
	traveltime_llr = cacheShortestTravelTimes(numlocs, prearcs, "llr time")
	orderArcSet, orderArcSet_space, A_plus_i, A_minus_i = Dict(), Dict(), Dict(), Dict()

	for i in orders
		orderArcSet[i] = [dummyarc]
		orderArcSet_space[i] = []
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
		for n2 in N_end
			arc_ext = extendedarcs[n2, extendednodes[destinationlocation, dummyendtime]]
			push!(orderArcSet[i], arc_ext)
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
			if traveltimefordelay_flag == 0
				t3 = traveltime_rdd[arcdest, endloc]
			elseif traveltimefordelay_flag >= 1
				t3 = traveltime_llr[arcdest, endloc]
			end

			#May want to toggle between these
			if ktype_flag == "abs"
				deadline = starttime + shortesttriptimes[i] + k
			elseif ktype_flag == "pct"
				deadline = starttime + shortesttriptimes[i] + k * shortesttriptimes[i]
			elseif ktype_flag == "min24"
				deadline = max(starttime + shortesttriptimes[i] + k * shortesttriptimes[i], starttime + 24)
			end

			for t in 0:tstep:horizon-t2
				if (starttime + t1 <= t) & (t + t2 + t3 <= deadline)
					push!(orderArcSet[i], arcs[nodes[arcorig, t], nodes[arcdest, t + t2]])
					if arcorig != arcdest
						push!(orderArcSet_space[i], arcs[nodes[arcorig, t], nodes[arcdest, t + t2]])
					end				
				end
			end
		end

	end

	#Create A_plus and A_minus lists
	for i in orders, n in 1:numnodes, a in A_plus[n]
		if (a in orderArcSet[i]) & !(a in A_plus_i[i,n])
			push!(A_plus_i[i,n], a)
		end
	end
	for i in orders, n in 1:numnodes, a in A_minus[n]
		if (a in orderArcSet[i]) & !(a in A_minus_i[i,n])
			push!(A_minus_i[i,n], a)
		end
	end

	return orderArcSet, orderArcSet_space, A_plus_i, A_minus_i

end

#---------------------------------------------------------------------------------------#

function initializeorderarcsets(k, orders, originloc, destloc, Origin, Destination, shortesttriptimes)

	if k == -1
		orderarcset, orderarcset_space, A_minus_i, A_plus_i = orderarcs_dummyonly(orders, originloc, destloc, Origin, Destination)
	else
		orderarcset, orderarcset_space, A_minus_i, A_plus_i = orderarcs_warmstart(k, orders, originloc, destloc, Origin, Destination, shortesttriptimes)
	end

	magarcs = (A=orderarcset, A_space=orderarcset_space, A_minus=A_minus_i, A_plus=A_plus_i, available=[])

	return magarcs

end