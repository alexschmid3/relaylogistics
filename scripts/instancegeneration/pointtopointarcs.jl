
function pointotpointarcs(orders, Origin, Destination)
	
	#For now, we are using shortest path regardless of driver availability (this would be more time consuming to do)
	traveltime_rdd = cacheShortestTravelTimes(numlocs, prearcs, "rdd time")
	#traveltime_llr = cacheShortestTravelTimes(numlocs, prearcs, "llr time")
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

    #Add all possible order delivery arcs
    for i in orders, t in nodesLookup[Origin[i][1]][2]:tstep:horizon
        deliverytime = arcLength[originloc[i], destloc[i]] #traveltime_rdd[originloc[i], destloc[i]] + (24-shiftlength)*floor(traveltime_rdd[originloc[i], destloc[i]] / shiftlength - 1e-4)
        arcorig, arcdest, arctraveltime = originloc[i], destloc[i], deliverytime
        arcfinishtime = t + arctraveltime > horizon ? dummyendtime : t + arctraveltime 
		
		push!(orderarcset[i], extendedarcs[extendednodes[arcorig, t], extendednodes[arcdest, arcfinishtime]])
        push!(orderarcset_space[i], extendedarcs[extendednodes[arcorig, t], extendednodes[arcdest, arcfinishtime]])
    end
	
	#Create A_plus and A_minus lists
	for i in orders, n in 1:extendednumnodes, a in A_plus[n]
		if (a in orderarcset[i]) & !(a in A_plus_i[i,n])
			push!(A_plus_i[i,n], a)
		end
	end
	for i in orders, n in 1:extendednumnodes, a in A_minus[n]
		if (a in orderarcset[i]) & !(a in A_minus_i[i,n])
			push!(A_minus_i[i,n], a)
		end
	end
	
	#Sort arcs by time (helpful for accelerating shortest path algorithms)
    for i in orders
		orderarcset[i] = union(sort(setdiff(orderarcset[i], dummyarc), by = x -> nodesLookup[arcLookup[x][1]][2]), dummyarc)
	end
	
	return orderarcset, orderarcset_space, A_plus_i, A_minus_i

end

#---------------------------------------------------------------------------------------#

