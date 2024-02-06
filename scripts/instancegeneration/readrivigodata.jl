using CSV, Random, DataFrames, Dates

#Read in location coordinates
function readlocations(filename, numlocs)

	data = CSV.read(filename, DataFrame)

	hubs_x = []
	hubs_y = []

	numlocs = min(maxlocs, size(data)[1])
	
	hubsLookup, hubsReverseLookup, hubsTravelTimeIndex = Dict(), Dict(), Dict()
	for i in 1:size(data)[1]
		push!(hubs_x, data[!,5][i])
		push!(hubs_y, data[!,6][i])
		hubsLookup[data[!,8][i]] = data[!,2][i]
		hubsReverseLookup[data[!,2][i]] = data[!,8][i]
		hubsTravelTimeIndex[data[!,1][i]] = data[!,8][i]
	end
	
	hubCoords = hcat(hubs_x, hubs_y)
	
	return hubCoords, hubsLookup, hubsReverseLookup, hubsTravelTimeIndex, numlocs

end

#---------------------------------------------------------------------------------------#

#Create time-space network of nodes
function timespacentwk(locationcount, timestep, timehorizon)
	nodes = Dict()
	nodesLookup = Dict()
	N_0 = []
	N_end = []
	index = 1
	for t in 0:timehorizon/timestep
		for l in 1:locationcount
			nodes[l,t*timestep] = index
			nodesLookup[index] = (l,t*timestep)
			if t == 0
				push!(N_0, index)
			elseif t == timehorizon/timestep
				push!(N_end, index)
			end
			index += 1
		end	
	end
	
	numnodes = length(nodes)
	
	return nodes, nodesLookup, N_0, N_end, numnodes
end

#---------------------------------------------------------------------------------------#

#Uniformly distribute trucks across the locations 
function truckdistribution(truckcount, locationcount, N_0, N_end)

	m_0 = []

	for n in N_0
		push!(m_0, floor(truckcount/locationcount))
	end 
	numtrucks_left = truckcount - floor(truckcount/locationcount)*locationcount
	for i in randperm(locationcount)[1:convert(Int64,numtrucks_left)]
		m_0[i] += 1
	end
	
	#Create ending truck counts
	m_end = Dict()
	for n in N_end
		m_end[n] = floor(0.9*truckcount/locationcount)
	end

	return m_0, m_end
	
end

#---------------------------------------------------------------------------------------#

function readarcs(filename, hubdistfilename, tstep, numlocs, hubsTravelTimeIndex, roundup_flag, excludeoutliers_flag, hubsReverseLookup, googlemapstraveltimes_flag)

	data = CSV.read(filename, DataFrame)

	#Find distances between pairs of nodes
	hubdistdata = CSV.read(hubdistfilename, DataFrame)

	distances = Dict()
	for i in 1:size(hubdistdata)[1]
		orig, dest = hubsReverseLookup[hubdistdata[i,1]], hubsReverseLookup[hubdistdata[i,2]]
		if (orig <= numlocs) & (dest <= numlocs)
			distances[orig, dest] = hubdistdata[i,3]
		end
	end

	#Create "pre-arcs": (loc1, loc2, travel_time) 
	prearcs = []
	arcLength, arcLength_raw = Dict(), Dict()

	if roundup_flag == 1
		for row in 1:size(data)[1]
			origin = hubsTravelTimeIndex[data[row,1]]
			destination = hubsTravelTimeIndex[data[row,2]]
			syntheticsymmetricarc_flag = data[row,7]

			#Add to pre-arcs if both origin and destination are in the list of locs for this instance
			if (origin <= numlocs) & (destination <= numlocs) & (syntheticsymmetricarc_flag <= includesymmetricarcs_flag)

				#Choose the travel time corresponding to the flags that are set (include outliers and/or use google maps travel times)
				if (googlemapstraveltimes_flag == 0) & (excludeoutliers_flag == 0)
					rawtt = data[row,5]
				elseif (googlemapstraveltimes_flag == 0) & (excludeoutliers_flag == 1)
					rawtt = data[row,6]
				elseif (googlemapstraveltimes_flag == 1) & (excludeoutliers_flag == 0)
					rawtt = max(data[row,5], data[row,8])
				elseif (googlemapstraveltimes_flag == 1) & (excludeoutliers_flag == 1)
					rawtt = max(data[row,6], data[row,8])
				elseif (googlemapstraveltimes_flag == 1) & (excludeoutliers_flag == 1) & (ensureconnectivity_flag == 1)
					if (shiftlength < data[row,6] < shiftlength + tstep) & (data[row,8] < shiftlength)
						rawtt = (data[row,6] + data[row,8])/2
					else
						rawtt = max(data[row,6], data[row,8])
					end
				end

				#Add travel time to list of pre-arcs and arcLength lists
				push!(prearcs,(origin, destination, ceil(rawtt/tstep) * tstep, rawtt, distances[origin, destination] ))
				arcLength[origin, destination] = ceil(rawtt/tstep) * tstep
				arcLength_raw[origin, destination] = rawtt
			end
		end
	elseif roundup_flag == 0
		for row in 1:size(data)[1]
			origin = hubsTravelTimeIndex[data[row,1]]
			destination = hubsTravelTimeIndex[data[row,2]]
			syntheticsymmetricarc_flag = data[row,7]
			
			#Add to pre-arcs if both origin and destination are in the list of locs for this instance
			if (origin <= numlocs) & (destination <= numlocs) & (syntheticsymmetricarc_flag <= includesymmetricarcs_flag)
				#Choose the travel time corresponding to the flags that are set (include outliers and/or use google maps travel times)
				if (googlemapstraveltimes_flag == 0) & (excludeoutliers_flag == 0)
					rawtt = data[row,5]
				elseif (googlemapstraveltimes_flag == 0) & (excludeoutliers_flag == 1)
					rawtt = data[row,6]
				elseif (googlemapstraveltimes_flag == 1) & (excludeoutliers_flag == 0)
					rawtt = max(data[row,5], data[row,8])
				elseif (googlemapstraveltimes_flag == 1) & (excludeoutliers_flag == 1)
					rawtt = max(data[row,6], data[row,8])
				end

				#Add travel time to list of pre-arcs and arcLength lists
				push!(prearcs,(origin, destination, floor(rawtt/tstep) * tstep, rawtt, distances[origin, destination] ))
				arcLength[origin, destination] = floor(rawtt/tstep) * tstep
				arcLength_raw[origin, destination] = rawtt
			end
		end
	end
	
	return prearcs, arcLength, arcLength_raw

end

#---------------------------------------------------------------------------------------#

function arccreation(prearcs, horizon, tstep, numnodes, nodes, numlocs)

	#Initialize output sets/dictionaries
	arcs = Dict()
	arcLookup = Dict()
	A_plus = Dict()
	A_minus = Dict()
	A_space = []
	A_plus_time = Dict()
	A_minus_time = Dict()
	A_minus_space = Dict()  
	A_plus_space = Dict()    
	index = 1
	truetraveltime = Dict()
	
	for node in 1:numnodes
		A_plus[node] = []
		A_minus[node] = []
		A_plus_time[node] = []
		A_minus_time[node] = []
		A_minus_space[node] = []
		A_plus_space[node] = []
	end
	
	#Create arcs
	for arc in prearcs
		traveltime = arc[3]
		for t in 0:tstep:horizon - traveltime
			startnode = nodes[arc[1],t]
			endnode = nodes[arc[2],t+arc[3]]
		
			arcs[(startnode,endnode)] = index
			arcLookup[index] = (startnode,endnode)
			
			push!(A_plus[startnode], index)
			push!(A_minus[endnode], index)
			push!(A_space, index)
			push!(A_minus_space[endnode], index)
			push!(A_plus_space[startnode], index)

			truetraveltime[index] = arc[4]
				
			index += 1
		end
	end
	
	for loc in 1:numlocs, t in 0:tstep:horizon - tstep
		startnode = nodes[loc,t]
		endnode = nodes[loc,t + tstep]

		arcs[(startnode,endnode)] = index
		arcLookup[index] = (startnode,endnode)
		
		push!(A_plus[startnode], index)
		push!(A_minus[endnode], index)	
		push!(A_plus_time[startnode], index)
		push!(A_minus_time[endnode], index)

		truetraveltime[index] = tstep
		
		index += 1
	end
	
	arccount = length(arcs)

	return arcs, arcLookup, A_plus, A_minus, A_space, A_plus_time, A_minus_time, A_minus_space, A_plus_space, arccount, truetraveltime
	
end

#---------------------------------------------------------------------------------------#

function prearcPreproc(prearcs, numlocs)

	P_plus = Dict()
	for l in 1:numlocs
		P_plus[l] = []
	end

	for pa in prearcs, l in 1:numlocs
		if pa[1] == l
			push!(P_plus[l], pa)
		end
	end
	
	return P_plus
	
end

#---------------------------------------------------------------------------------------#

function findShortestPath(loc1, loc2, numlocs, prearcs, P_plus, arcweighttype)

	#Initialize shortest path algorithm (Dijkstra's)
	visitednodes = zeros(numlocs)
	currdistance = repeat([999999.0],outer=[numlocs])
	currdistance[loc1] = 0
	currloc = loc1
	nopathexistsflag = 0
	prevloc = zeros(numlocs)

	#Find shortest path from loc1 to loc2
	while (visitednodes[loc2] == 0) & (nopathexistsflag == 0)

		#Assess all neighbors of current node
		for arc in P_plus[currloc]
			l = arc[2]
			if visitednodes[l] == 0
				if arcweighttype == "rdd time"
					newdist = currdistance[currloc] + arc[3]
				elseif (arcweighttype == "llr time") & (l != loc2)
					newdist = currdistance[currloc] + arc[3]
				elseif (arcweighttype == "llr time") & (l == loc2)
					newdist = currdistance[currloc] + arc[4]
				elseif arcweighttype == "raw time"
					newdist = currdistance[currloc] + arc[4]
				elseif arcweighttype == "dist"
					newdist = currdistance[currloc] + arc[5]
				end
			
				if newdist < currdistance[l]
					currdistance[l] = newdist
					prevloc[l] = currloc
				end
			end
		end

		#Mark the current node as visited
		visitednodes[currloc] = 1

		#Update the current node 
		currdistance_unvisited = deepcopy(currdistance)
		for l in 1:numlocs
			if visitednodes[l] == 1
				currdistance_unvisited[l] = 999999
			end
		end
		currloc = argmin(currdistance_unvisited)

		#If all remaining nodes have tentative distance of 999999 and the algorithm has not terminated, then there is no path from origin to destination
		if minimum(currdistance_unvisited) == 999999
			nopathexistsflag = 1
		end

	end

	#Format the shortest path output
	shortestpatharcs_rev = []
	node = loc2
	while node != loc1
		push!(shortestpatharcs_rev, (Int(prevloc[node]), node))
		node = Int(prevloc[node])
	end
	shortestpatharcs = reverse(shortestpatharcs_rev) 

	return currdistance[loc2], shortestpatharcs

end

#---------------------------------------------------------------------------------------#

function cacheShortestTravelTimes(numlocs, prearcs, arcweighttype)
	
	P_plus = prearcPreproc(prearcs, numlocs)
	shortestTravelTime = Dict()

	for loc1 in 1:numlocs, loc2 in 1:numlocs

		if loc1 == loc2
			shortestTravelTime[loc1, loc2] = 0
		else
			tt, sparcs = findShortestPath(loc1, loc2, numlocs, prearcs, P_plus, arcweighttype)
			shortestTravelTime[loc1, loc2] = tt
		end

	end

	return shortestTravelTime

end

#---------------------------------------------------------------------------------------#

function cacheShortestDistance(numlocs, prearcs)
	
	P_plus = prearcPreproc(prearcs, numlocs)
	shortestPathArcs = Dict()
	shortestDistance = Dict()

	for loc1 in 1:numlocs, loc2 in 1:numlocs

		if loc1 == loc2
			shortestDistance[loc1, loc2] = 0
		else
			#shortestDistance[loc1, loc2] = findShortestPath(loc1, loc2, numlocs, prearcs, P_plus, "dist")
			tt, arclist = findShortestPath(loc1, loc2, numlocs, prearcs, P_plus, "dist")
			shortestDistance[loc1, loc2] = tt
			shortestPathArcs[loc1, loc2] = arclist
		end

	end

	return shortestDistance, shortestPathArcs

end

#---------------------------------------------------------------------------------------#

function generateorderlist(lh_filename, vnt_filename, ordermaxcap, numlocs)

	data_agg = CSV.read(lh_filename, DataFrame)

	orderDictHelper = Dict()
	for i in 1:size(data_agg)[1]
		orig, dest = data_agg[!,26][i], data_agg[!,27][i]
		psseq_raw = data_agg[i,8]
		psseq = split(psseq_raw, "-")
		orderid = data_agg[!,1][i]

		pickup_ts = DateTime(1970) + Dates.Millisecond(data_agg[!,6][i])
		deliv_ts = DateTime(1970) + Dates.Millisecond(data_agg[!,5][i])

		#start_avail_ts = floor(pickup_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8)			
		start_avail_ts = floor(pickup_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8) + 6*floor(Dates.Millisecond(floor(Dates.value(pickup_ts - (floor(pickup_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8)))/6)), Dates.Hour)
		end_avail_ts = start_avail_ts + Dates.Day(1)
		start_due_ts = floor(deliv_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8)
		end_due_ts = start_due_ts + Dates.Day(1) 

		#Check whether all intermediate nodes from the Rivigo pitstop sequence are included in the subset of locs
		intermedlocs_flag = 0
		hubsList = collect(values(hubsLookup))
		for ps in psseq
			if ps in hubsList
				loc = hubsReverseLookup[ps]
				if loc > numlocs
					intermedlocs_flag = 1
					break
				end	
			else
				intermedlocs_flag = 1
				break
			end
		end

		if (orig != dest) & (1 <= orig <= numlocs) & (1 <= dest <= numlocs) & (intermedlocs_flag == 0) #& (orderwindowstart <= start_avail_ts <= orderwindowend) 
			try 
				push!(orderDictHelper[start_avail_ts], orderid)
			catch 
				orderDictHelper[start_avail_ts] = [orderid]
			end
		end
	end

	#timeblocklist = sort(collect(keys(orderDictHelper)))
	timeblocklist = collect(keys(orderDictHelper))

	#Randomly permute the orders in each time block
	for block in timeblocklist
		orderDictHelper[block] = orderDictHelper[block][randperm(length(orderDictHelper[block]))]
	end

	includelist = []
	weeklycounterdict = Dict()
	for j in 1:ordermaxcap
		for block in timeblocklist
			start_avail_wk = floor(block - Dates.Hour(8), Dates.Week) + Dates.Hour(8) 
			if orderDictHelper[block] != []
				chosenorder = popfirst!(orderDictHelper[block])
				push!(includelist, chosenorder)
				try 
					push!(weeklycounterdict[start_avail_wk], chosenorder)
				catch 
					weeklycounterdict[start_avail_wk] = [chosenorder]
				end
			end
		end
	end

	return includelist

end

#---------------------------------------------------------------------------------------#

function pullorders_initrivigoroutes(lh_filename, vnt_filename, maxorders, orderwindowstart, orderwindowend, tstep, horizon, prearcs, numlocs, timedelta, includelist)

	traveltime = cacheShortestTravelTimes(numlocs, prearcs, "rdd time")

	data_agg = CSV.read(lh_filename, DataFrame)
	data_byleg = CSV.read(vnt_filename, DataFrame)
	
	originloc, destloc, available, duedate, orderidlist, psseqlist, trueoriginloc, orderinprog = [], [], [], [], [], [], [], []

	currentLocation = Dict()
	for id in 0:maximum(data_byleg[!,5])
		currentLocation[id] = (99999, DateTime(2099, 1, 1, 8))
		#currentLocation[id] = (99999, 99999)
	end

	#Find the current location for any orders that happen within our time period of interest
	for leg in 1:size(data_byleg)[1]
		orderid = data_byleg[!,5][leg]
		hubcode = data_byleg[!,2][leg]
		if !ismissing(data_byleg[!,3][leg]) & !(hubcode == "BLSP1") & !(hubcode == "BWLP1")
			in_ts = DateTime(1970) + Dates.Millisecond(data_byleg[!,3][leg])
			if (currentLocation[orderid] == (99999, DateTime(2099, 1, 1, 8))) & (in_ts >= orderwindowstart) & (in_ts <= orderwindowend + Dates.Hour(timedelta))
				#starttime = tstep* ceil(((in_ts - orderwindowstart) / (Millisecond(1) * 1000 * 3600)) / tstep)
				currentLocation[orderid] = (hubsReverseLookup[hubcode], in_ts)
			end
		end
	end

	#Filter to find the orders 
	for i in 1:size(data_agg)[1]
		orig, dest = data_agg[!,26][i], data_agg[!,27][i]
		psseq_raw = data_agg[i,8]
		psseq = split(psseq_raw, "-")
		orderid = data_agg[!,1][i]
		if orderid in includelist
			
			pickup_ts = DateTime(1970) + Dates.Millisecond(data_agg[!,6][i])
			deliv_ts = DateTime(1970) + Dates.Millisecond(data_agg[!,5][i])

			#start_avail_ts = floor(pickup_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8)			
			start_avail_ts = floor(pickup_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8) + timedelta*floor(Dates.Millisecond(floor(Dates.value(pickup_ts - (floor(pickup_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8)))/timedelta)), Dates.Hour)
			end_avail_ts = start_avail_ts + Dates.Day(1)
			start_due_ts = floor(deliv_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8)
			end_due_ts = start_due_ts + Dates.Day(1) 

			#Check whether all intermediate nodes from the Rivigo pitstop sequence are included in the subset of locs
			intermedlocs_flag = 0
			hubsList = collect(values(hubsLookup))
			for ps in psseq
				if ps in hubsList
					loc = hubsReverseLookup[ps]
					if loc > numlocs
						intermedlocs_flag = 1
						break
					end	
				else
					intermedlocs_flag = 1
					break
				end
			end

			if (orig != dest) & (1 <= orig <= numlocs) & (1 <= dest <= numlocs) & (intermedlocs_flag == 0) & (orderwindowstart <= start_avail_ts <= orderwindowend) 
			#if (orig != dest) & (1 <= orig <= numlocs) & (1 <= dest <= numlocs) & (orderwindowstart <= start_avail_ts <= orderwindowend) & (orderwindowstart <= end_due_ts <= orderwindowend ) & (start_avail_ts <= start_due_ts)

				start_avail = (start_avail_ts - orderwindowstart) / (Millisecond(1) * 1000 * 3600)
				#end_avail = (end_avail_ts - orderwindowstart) / (Millisecond(1) * 1000 * 3600)
				end_avail = horizon
				#start_due = (start_due_ts - orderwindowstart) / (Millisecond(1) * 1000 * 3600)
				start_due = min(horizon, start_avail + traveltime[orig, dest])
				#end_due = (end_due_ts - orderwindowstart) / (Millisecond(1) * 1000 * 3600)
				end_due = horizon

				push!(originloc, orig)
				push!(destloc, dest)
				push!(available, (start_avail, end_avail))
				push!(duedate, (start_due, end_due))
				push!(orderidlist, orderid)
				push!(psseqlist, psseq)
				push!(trueoriginloc, orig)

			elseif (currentLocation[orderid][1] != orig) & (orig != dest) & (currentLocation[orderid][1] != dest) & (1 <= orig <= numlocs) & (1 <= dest <= numlocs) & (intermedlocs_flag == 0) & (1 <= currentLocation[orderid][1] <= numlocs) & (orderwindowstart <= currentLocation[orderid][2] <= orderwindowend) 
			#elseif (orig != dest) & (1 <= orig <= numlocs) & (1 <= dest <= numlocs) & (1 <= currentLocation[orderid][1] <= numlocs) & (orderwindowstart <= currentLocation[orderid][2] <= orderwindowend) & (orderwindowstart <= end_due_ts <= orderwindowend) & (start_avail_ts <= start_due_ts)
			
				start_avail_ts2 = floor(currentLocation[orderid][2] - Dates.Hour(8), Dates.Day) + Dates.Hour(8) + timedelta*floor(Dates.Millisecond(floor(Dates.value(currentLocation[orderid][2] - (floor(currentLocation[orderid][2] - Dates.Hour(8), Dates.Day) + Dates.Hour(8)))/timedelta)), Dates.Hour)
				start_avail = (start_avail_ts2 - orderwindowstart) / (Millisecond(1) * 1000 * 3600)
				#end_avail = (start_avail_ts2 + Dates.Day(1) - orderwindowstart) / (Millisecond(1) * 1000 * 3600)
				#end_avail = horizon
				end_avail = start_avail

				start_due = min(horizon, start_avail + traveltime[currentLocation[orderid][1], dest])
				end_due = horizon

				push!(originloc, currentLocation[orderid][1])
				push!(destloc, dest)
				push!(available, (start_avail, end_avail))
				push!(duedate, (start_due, end_due))
				push!(orderidlist, orderid)
				push!(psseqlist, psseq)
				push!(trueoriginloc, orig)
				push!(orderinprog, length(originloc))

			end
		end
	end

	numorders = min(maxorders, length(originloc))

	return numorders, originloc[1:numorders], destloc[1:numorders], available[1:numorders], duedate[1:numorders], orderidlist[1:numorders], psseqlist[1:numorders], trueoriginloc[1:numorders], [i for i in orderinprog if i <= numorders]

end

#---------------------------------------------------------------------------------------#

function pullorders_rivigoroutes(lh_filename, vnt_filename, maxorders, orderwindowstart, orderwindowend, currentdatetime, tstep, horizon, prearcs, numlocs, timedelta, includelist)
        
    #excludelist = [5625,5703,6207,6280,6325,6330,6338,6361,6385,6400,6412,6428,6429,6447,6448,6455,6489,6508,6553,6559,6577,6601,6613,6632,6634,6651,6664,6686,6690,6714,6727,6845,6879,6940,6948,7035,7066,7073,7074,7097,7098,7102,7113,7140,7154,7164,7194,7264,7268,7287,7295,7310,7337,7338,7345,7364,7438,7454,7463,7529,7544,7557,7571,7601,7603,7636,7639,7653,7662,7676,7679,7689,7736,7741,7793,7818,7824,7866,7885,7921,7976,7980,8035,8065,8078,8119,8121,8145,8177,8180,8189,8372,8389,8406,8414,8456,8467,8483,8495,8516,8611,8640,8644,8666,8701,8707,8728,8742,8770,8786,8800,8815,8823,8824,8836,8852,8916,8937,8942,8966,8997,8998,9002,9010,9011,9025,9043,9056,9069,9096,9140,9164,9219,9237,9254,9262,9267,9355,9374,9400,9424,9428,9447,9464,9480,9504,9516,9544,9546,9558,9566,9567,9574,9647,9683,9689,9701,9734,9743,9744,9778,9805,9840,9863,9879,9893,9905,9924,9967,9976,9992,9997]
	#P_plus = prearcPreproc(prearcs, numlocs)
	traveltime = cacheShortestTravelTimes(numlocs, prearcs, "rdd time")

	data_agg = CSV.read(lh_filename, DataFrame)
	data_byleg = CSV.read(vnt_filename, DataFrame)
	
	originloc, destloc, available, duedate, orderidlist, psseqlist, trueoriginloc, orderinprog = [], [], [], [], [], [], [], []

	currentLocation = Dict()
	for id in 0:maximum(data_byleg[!,5])
		currentLocation[id] = (99999, DateTime(2099, 1, 1, 8))
		#currentLocation[id] = (99999, 99999)
	end

	#Find the current location for any orders that happen within our time period of interest
	for leg in 1:size(data_byleg)[1]
		orderid = data_byleg[!,5][leg]
		hubcode = data_byleg[!,2][leg]
		if !ismissing(data_byleg[!,3][leg]) & !(hubcode == "BLSP1") & !(hubcode == "BWLP1")
			in_ts = DateTime(1970) + Dates.Millisecond(data_byleg[!,3][leg])
			if (currentLocation[orderid] == (99999, DateTime(2099, 1, 1, 8))) & (in_ts >= orderwindowstart) & (in_ts <= orderwindowend + Dates.Hour(timedelta))
				#starttime = tstep* ceil(((in_ts - orderwindowstart) / (Millisecond(1) * 1000 * 3600)) / tstep)
				currentLocation[orderid] = (hubsReverseLookup[hubcode], in_ts)
			end
		end
	end

	#Filter to find the orders 
	for i in 1:size(data_agg)[1]
		orig, dest = data_agg[i,26], data_agg[i,27]
		psseq_raw = data_agg[i,8]
		psseq = split(psseq_raw, "-")
		orderid = data_agg[i,1]
		if orderid in includelist
			pickup_ts = DateTime(1970) + Dates.Millisecond(data_agg[!,6][i])
			deliv_ts = DateTime(1970) + Dates.Millisecond(data_agg[!,5][i])

			#start_avail_ts = floor(pickup_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8)			
			start_avail_ts = floor(pickup_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8) + timedelta*floor(Dates.Millisecond(floor(Dates.value(pickup_ts - (floor(pickup_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8)))/timedelta)), Dates.Hour)
			end_avail_ts = start_avail_ts + Dates.Day(1)
			start_due_ts = floor(deliv_ts - Dates.Hour(8), Dates.Day) + Dates.Hour(8)
			end_due_ts = start_due_ts + Dates.Day(1) 

			#Check whether all intermediate nodes from the Rivigo pitstop sequence are included in the subset of locs
			intermedlocs_flag = 0
			hubsList = collect(values(hubsLookup))
			for ps in psseq
				if ps in hubsList
					loc = hubsReverseLookup[ps]
					if loc > numlocs
						intermedlocs_flag = 1
						break
					end	
				else
					intermedlocs_flag = 1
					break
				end
			end

			if (orig != dest) & (1 <= orig <= numlocs) & (1 <= dest <= numlocs) & (intermedlocs_flag == 0) & (orderwindowstart <= start_avail_ts <= orderwindowend) 
			#if (orig != dest) & (1 <= orig <= numlocs) & (1 <= dest <= numlocs) & (orderwindowstart <= start_avail_ts <= orderwindowend) & (orderwindowstart <= end_due_ts <= orderwindowend ) & (start_avail_ts <= start_due_ts)

				start_avail = (start_avail_ts - currentdatetime) / (Millisecond(1) * 1000 * 3600)
				#end_avail = (end_avail_ts - currentdatetime) / (Millisecond(1) * 1000 * 3600)
				end_avail = horizon
				#start_due = (start_due_ts - currentdatetime) / (Millisecond(1) * 1000 * 3600)
				start_due = min(horizon, start_avail + traveltime[orig, dest])
				#end_due = (end_due_ts - currentdatetime) / (Millisecond(1) * 1000 * 3600)
				end_due = horizon

				push!(originloc, orig)
				push!(destloc, dest)
				push!(available, (start_avail, end_avail))
				push!(duedate, (start_due, end_due))
				push!(orderidlist, orderid)
				push!(psseqlist, psseq)
				push!(trueoriginloc, orig)

			end
		end
	end

	numorders = min(maxorders, length(originloc))

	return numorders, originloc[1:numorders], destloc[1:numorders], available[1:numorders], duedate[1:numorders], orderidlist[1:numorders], psseqlist[1:numorders], trueoriginloc[1:numorders], [i for i in orderinprog if i <= numorders]

end

#---------------------------------------------------------------------------------------#

function formatorders(numorders, originloc, destloc, available, duedate, tstep)

	Origin, Destination = Dict(), Dict()

	for i in 1:numorders
		Origin[i] = []
		Destination[i] = []

		for t in available[i][1]:tstep:available[i][2]
			push!(Origin[i], nodes[originloc[i], t])
		end

		for t in duedate[i][1]:tstep:duedate[i][2]
			push!(Destination[i], nodes[destloc[i], t])
		end
	end

	return Origin, Destination

end

#---------------------------------------------------------------------------------------#

function readdrivers(filename, maxdrivers, numlocs, nodes, horizon)

	data = CSV.read(filename, DataFrame)
	
	driverStartNodes, driverEndNodes, driverHomeLocs, assignedDrivers = Dict(), Dict(), Dict(), Dict()
	realdrivers, adjusteddrivers = Dict(), Dict()
	for l in 1:numlocs
		assignedDrivers[l] = []
	end

	totaldrivers = 0
	for i in 1:size(data)[1]
		loc = data[i,2]
		drivers = data[i,3]
		if loc <= numlocs
			totaldrivers += drivers
			realdrivers[loc] = drivers
		end
	end

	adjtotal = 0
	if totaldrivers >= maxdrivers
		for loc in 1:numlocs
			adjusteddrivers[loc] = max(1, floor(maxdrivers * realdrivers[loc] / totaldrivers))
			adjtotal += max(1, floor(maxdrivers * realdrivers[loc] / totaldrivers))
		end

		remaining = maxdrivers - adjtotal
		if remaining < 0 
			println("You allocated too many drivers :(")
		end

		for loc in randperm(numlocs)[1:convert(Int64,remaining)]
			adjusteddrivers[loc] += 1
		end
	else
		adjusteddrivers = realdrivers
	end

	driver = 1
	for loc in 1:numlocs
		for d in 1:adjusteddrivers[loc]
			driverHomeLocs[driver] = loc
			driverStartNodes[driver] = nodes[loc, 0]
			driverEndNodes[driver] = nodes[loc, horizon]
			push!(assignedDrivers[loc], driver)

			driver += 1
		end
	end

	driverindices = 1:convert(Int64,sum(adjusteddrivers[loc] for loc in 1:numlocs))

	return driverindices, driverStartNodes, driverEndNodes, driverHomeLocs, assignedDrivers

end

#---------------------------------------------------------------------------------------#

function flowbalancenodesets_i(orderindices, nodecount, OrigNodes, DestNodes)
	N_flow_i = Dict()
	for i in orderindices
		N_flow_i[i] = []
	end

	for i in orderindices
		for n in 1:nodecount
			if n in OrigNodes[i]
				1+1
			elseif n in DestNodes[i]
				1+1
			else
				push!(N_flow_i[i], n)
			end
		end
	end
	
	return N_flow_i

end
	
#---------------------------------------------------------------------------------------#

function flowbalancenodesets_td(driverindices, nodecount, driverStartLocations, N_start, N_finish)

	N_flow_t = []
	for n in 1:nodecount
		if n in N_start	
			1+1
		elseif n in N_finish
			1+1
		else 
			push!(N_flow_t, n)
		end
	end

	N_flow_d = Dict()
	for d in driverindices
		N_flow_d[d] = []
	end

	for d in driverindices
		for n in 1:nodecount
			if n in driverStartLocations[d]	
				1+1
			elseif n in N_finish
				1+1
			else 
				push!(N_flow_d[d], n)
			end
		end
	end

	return N_flow_t, N_flow_d

end

#---------------------------------------------------------------------------------------#

function driverArcSetsByDriver_overnight(numlocs, numarcs, numnodes, prearcs, drivers, tstep, horizon, nodes, arcs, assignedDrivers, A_minus, A_plus, T_off, drivershift, driverHomeLocs, T_off_0, shiftlength)

	#P_plus = prearcPreproc(prearcs, numlocs)
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

function driverArcSetsByDriver_nonrelay(numlocs, numarcs, numnodes, prearcs, drivers, tstep, horizon, nodes, arcs, assignedDrivers, A_minus, A_plus, T_off, drivershift, driverHomeLocs, T_off_0, shiftlength)

	#P_plus = prearcPreproc(prearcs, numlocs)
	traveltime = cacheShortestTravelTimes(numlocs, prearcs, "rdd time")
	homeArcSet, homeArcSet_space, availableDrivers, A_plus_d, A_minus_d = Dict(), Dict(), Dict(), Dict(), Dict()
	fullhourlist = [t for t in 0:tstep:horizon-tstep]

	for d in drivers
		homeArcSet[d] = []
		homeArcSet_space[d] = []
	end
	for a in 1:extendednumarcs
		availableDrivers[a] = []
	end
	for d in drivers, n in 1:numnodes
		A_plus_d[d,n] = []
		A_minus_d[d,n] = []
	end
	for d in drivers, n in numnodes+1:extendednumnodes
		A_minus_d[d,n] = []
	end

	prearcs_aug = deepcopy(prearcs)
	for l in 1:numlocs
		push!(prearcs_aug, (l, l, tstep, tstep))
	end

	closelocs = Dict()
	for d in drivers
		closelocs[d] = [l for l in 1:numlocs]
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

		#Add arc to the driver's arc list
		for t in setdiff([t4 for t4 in fullhourlist if t4 >= earliesttime], T_off[drivershift[d]])
			#Arc must finish before the end of the horizon and before the "next" off hour of the driver
			if t + t2 <= [t4 for t4 in union(T_off[drivershift[d]], horizon) if t4 > t][1]
				push!(homeArcSet[d], arcs[nodes[orig, t], nodes[dest, t + t2]])
				if orig != dest
					push!(homeArcSet_space[d], arcs[nodes[orig, t], nodes[dest, t + t2]])
				end				
				push!(availableDrivers[arcs[nodes[orig, t], nodes[dest, t + t2]]], d)
			end
		end
		
	end

	#Add stay-at-home arcs for driver's off hours
	for d in drivers, t in setdiff(T_off[drivershift[d]], [horizon]), l in closelocs[d]
		push!(homeArcSet[d], arcs[nodes[l, t], nodes[l, t + tstep]])
		push!(availableDrivers[arcs[nodes[l, t], nodes[l, t + tstep]]], d)
	end

	#Add extended arcs 
	for d in drivers, n1 in N_end, n2 in numnodes+1:extendednumnodes
		a = extendedarcs[n1,n2]
		push!(homeArcSet[d], a)
		push!(A_minus_d[d,n2], a)
		#push!(A_plus_d[d,n1], a)
		push!(availableDrivers[a], d)
	end

	#Create A_plus and A_minus lists
	for d in drivers, n in 1:numnodes, a in A_plus[n]
		if a in homeArcSet[d]
			push!(A_plus_d[d,n], a)
		end
	end
	for d in drivers, n in 1:numnodes, a in A_minus[n]
		if a in homeArcSet[d]
			push!(A_minus_d[d,n], a)
		end
	end

	return homeArcSet, homeArcSet_space, availableDrivers, A_plus_d, A_minus_d, closelocs

end

#---------------------------------------------------------------------------------------#

#Create objective costs
function calcobjectivecosts(arccount, arcLookupDict, nodesLookupDict,  locations)
	
	c = []

	for i in 1:arccount
		originloc = nodesLookupDict[arcLookupDict[i][1]][1]
		destloc = nodesLookupDict[arcLookupDict[i][2]][1]
		dist = sqrt((locations[originloc,1]-locations[destloc,1])^2 + (locations[originloc,2]-locations[destloc,2])^2)
		
		if dist == 0
			push!(c, 0)
		else 
			push!(c, dist)
		end 
	end
	
	return c

end

#---------------------------------------------------------------------------------------#

#Read objective costs using GreatCircleDistances
function readobjectivecosts(filename, numarcs, numlocs, hubsReverseLookup, nodesLookup, arcLookup)
	
	data = CSV.read(filename, DataFrame)

	distances = Dict()

	for i in 1:size(data)[1]
		orig, dest = hubsReverseLookup[data[i,1]], hubsReverseLookup[data[i,2]]
		if (orig <= numlocs) & (dest <= numlocs)
			distances[orig, dest] = data[i,3]
		end
	end

	c = [0.0 for a in 1:numarcs]

	for i in 1:numarcs
		oloc = nodesLookup[arcLookup[i][1]][1]
		dloc = nodesLookup[arcLookup[i][2]][1]
				
		c[i] += distances[oloc, dloc]

	end
	
	return c

end

#---------------------------------------------------------------------------------------#

#Cluster pitstops within a certain radius for DP
function clusterpitstops(prearcs, dp_pitstopclusterradius)
	
	clusteredpitstops = Dict()

	for l in 1:numlocs
		clusteredpitstops[l] = []
	end

	for arc in prearcs			
		oloc, dloc, drivedist = arc[1], arc[2], arc[5]
		if drivedist <= dp_pitstopclusterradius
			push!(clusteredpitstops[oloc],dloc)
		end
	end
	
	return clusteredpitstops 

end

#---------------------------------------------------------------------------------------#

function yarcreduction(numarcs, availableDrivers, A_space, numnodes, A_plus, A_minus)

	A_hasdriver = []
	yupperbound = []
	A_hasdriver_space = []
	for a in 1:numarcs
		ub = 99999
		if a in A_space
			ub = length(availableDrivers[a])
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

function orderarcreduction(prearcs, shortesttriptime)
	
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
		for n2 in setdiff(N_end, Destination[i])
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

			for t in 0:tstep:horizon-t2
				if starttime + t1 <= t
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

	for i in orders
		orderArcSet[i] = union(sort(setdiff(orderArcSet[i], dummyarc), by = x -> nodesLookup[arcLookup[x][1]][2]), dummyarc)
	end

	return orderArcSet, orderArcSet_space, A_plus_i, A_minus_i

end

#---------------------------------------------------------------------------------------#

function orderarcreduction_unik(k, ktype_flag, prearcs, shortesttriptimes)
	
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

function createdrivershifts(drivers, shiftlength, tstep, drivershifttstep, alltimeswithinview)

	#Driver shifts
	T_off_Monday8am = []	
	baseshift = []
	for day in -1:7, t in 0:tstep:24-tstep
		if 0 <= t < shiftlength
			push!(baseshift, day*24 + t)
		end
	end
	if driveroffdays_flag == 1
		offdaypairs = [(d1,mod(d1,7)+1) for d1 in 1:7]
	elseif driveroffdays_flag == 0
		offdaypairs = [(-10,-10)]
	end	
	for off in offdaypairs, s in 0:drivershifttstep:24-drivershifttstep
		shiftedshift = []
		for t in 0:tstep:168-tstep
			if (t-s in baseshift) #& (0 <= t <= 168-tstep)
				push!(shiftedshift, t)
			elseif convert(Int64, floor(t/24)) + 1 in off
				push!(shiftedshift, t)
			end
		end
		extendedshift = []
		for wk in 1:convert(Int64, ceil(alltimeswithinview/168))
			extendedshift = vcat(extendedshift, [168*(wk-1) for y in 1:length(shiftedshift)] + shiftedshift)
		end
		push!(T_off_Monday8am, extendedshift)
	end

	numshifts = length(T_off_Monday8am)

	T_off = []
	for shift in T_off_Monday8am
		newshift = []
		for hr in shift
			if hr <= horizon #- tstep
				push!(newshift, hr)
			end
		end
		push!(T_off, newshift)
	end

	drivershift = []
	for d in drivers
		push!(drivershift, rand(1:numshifts))
	end

	T_off_0, T_off_constr = Dict(), Dict()
	for d in drivers
		T_off_0[d], T_off_constr[d] = [], []
		for t in 0:tstep:horizon
			if (t in T_off[drivershift[d]]) & !(t-tstep in T_off[drivershift[d]])
				push!(T_off_0[d], t)
				if (t + 24 <= horizon) & (intersect([t3 for t3 in t:tstep:t+24],T_off_0[d]) != [])
					push!(T_off_constr[d], t)  
				end
			end
		end
	end

	#Online remediation - this can't be dependent on T_off_Monday8am
	#Create a new set, T_off_ext, which is T_off adjusted for new online time horizon but extended 24*maxnightsaway past the time horizon end
	T_on_0 = Dict()
	for d in drivers
		T_on_0[d] = []
		for t in 0:tstep:horizon+24*maxnightsaway
			if (t in setdiff([t2 for t2 in 0:tstep:horizon+24*maxnightsaway], T_off_Monday8am[drivershift[d]])) & !(t-tstep in setdiff([t2 for t2 in 0:tstep:horizon+24*maxnightsaway], T_off_Monday8am[drivershift[d]]))
				push!(T_on_0[d], t)
			end
		end
	end

	return T_off_Monday8am, T_off, drivershift, T_off_0, T_off_constr, numshifts, T_on_0

end

#---------------------------------------------------------------------------------------#

function createonlineshifts(shiftlength)

	#Driver shifts
	T_off_Monday8am_online = []	
	baseshift = []
	for day in -1:7, t in 0:tstep:24-tstep
		if 0 <= t < shiftlength
			push!(baseshift, day*24 + t)
		end
	end
	if driveroffdays_flag == 1
		offdaypairs = [(d1,mod(d1,7)+1) for d1 in 1:7]
	elseif driveroffdays_flag == 0
		offdaypairs = [(-10,-10)]
	end	
	for off in offdaypairs, s in 0:tstep:24-tstep
		shiftedshift = []
		for t in 0:tstep:168-tstep
			if (t-s in baseshift) #& (0 <= t <= 168-tstep)
				push!(shiftedshift, t)
			elseif convert(Int64, floor(t/24)) + 1 in off
				push!(shiftedshift, t)
			end
		end
		extendedshift = []
		for wk in 1:convert(Int64, ceil(alltimeswithinview/168))
			extendedshift = vcat(extendedshift, [168*(wk-1) for y in 1:length(shiftedshift)] + shiftedshift)
		end
		push!(T_off_Monday8am_online, extendedshift)
	end

	numonlineshifts = length(T_off_Monday8am_online)

	T_off_online = []
	onlineshiftlookup = Dict()
	index = 1
	for shift in T_off_Monday8am_online
		newshift = []
		for hr in shift
			if hr <= horizon #- tstep
				push!(newshift, hr)
			end
		end
		push!(T_off_online, newshift)
		onlineshiftlookup[newshift] = index
		index += 1
	end

	#Shift start times
	T_on_0 = Dict()
	for os in 1:numonlineshifts
		T_on = setdiff(0:tstep:horizon, T_off_online[os])
		T_on_0[os] = [t for t in T_on if !(t-tstep in T_on)]
	end
	
	return T_off_online, numonlineshifts, onlineshiftlookup, T_on_0

end

#---------------------------------------------------------------------------------------#

function createdrivershifts_nonrelay(shiftlength, tstep, drivershifttstep)

	#Driver shifts
	T_off_Monday8am = []	
	baseshift = []
	for day in -1:7, t in 0:tstep:24-tstep
		if 0 <= t < shiftlength
			push!(baseshift, day*24 + t)
		end
	end
	for s in 0:drivershifttstep:24-drivershifttstep
		shiftedshift = []
		for t in baseshift
			if 0 <= t+s <= 168-tstep
				push!(shiftedshift, t+s)
			end
		end
		extendedshift = []
		for wk in 1:convert(Int64, max(driverreturnhomeweeks, ceil(onlinetimehorizon/168) + 1, ceil(alltimeswithinview/168)))
			extendedshift = vcat(extendedshift, [168*(wk-1) for y in 1:length(shiftedshift)] + shiftedshift)
		end
		push!(T_off_Monday8am, extendedshift)
	end
	numshifts = length(T_off_Monday8am)

	T_off = []
	for shift in T_off_Monday8am
		newshift = []
		for hr in shift
			if hr <= horizon - tstep
				push!(newshift, hr)
			end
		end
		push!(T_off, newshift)
	end

	drivershift = []
	for d in drivers
		push!(drivershift, rand(1:numshifts))
	end

	T_off_0, T_off_constr = Dict(), Dict()
	for d in drivers
		T_off_0[d], T_off_constr[d] = [], []
		for t in 0:tstep:horizon
			if (t in T_off[drivershift[d]]) & !(t-tstep in T_off[drivershift[d]])
				push!(T_off_0[d], t)
				if (t + 24 <= horizon) & (intersect([t3 for t3 in t:tstep:t+24],T_off_0[d]) != [])
					push!(T_off_constr[d], t)  
				end
			end
		end
	end

	T_on_0 = Dict()
	for d in drivers
		T_on_0[d] = []
		for t in 0:tstep:horizon
			if (t in setdiff([t2 for t2 in 0:tstep:horizon], T_off[drivershift[d]])) & !(t-tstep in setdiff([t2 for t2 in 0:tstep:horizon], T_off[drivershift[d]]))
				push!(T_on_0[d], t)
			end
		end
	end


	return T_off_Monday8am, T_off, drivershift, T_off_0, T_off_constr, numshifts, T_on_0

end

#---------------------------------------------------------------------------------------#

function calcobjectivecosts(hubdistancesfilename)
	
	c = readobjectivecosts(hubdistancesfilename, numarcs, numlocs, hubsReverseLookup, nodesLookup, arcLookup)
	u = [0.0 for a in 1:numarcs]
	for a in 1:numarcs
		u[a] = taxicostpct*c[a]  
	end

	return c, u
end

#---------------------------------------------------------------------------------------#

println("Loaded Rivigo data helper functions")