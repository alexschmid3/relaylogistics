
include("convergenceplots.jl")

using Base.Iterators: partition

#---------------------------------------------------------------------------------------#

function initializearcsets(orders)
	
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

function initializearcsets_warmstart(k, ktype_flag)
	
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

#----------------------------------------------------------------------------------------#

function sparsemasterproblem(orderArcSet, orderArcSet_space, A_plus_i, A_minus_i, homeArcSet, A_minus_d, A_plus_d, availableDrivers)

	smp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(smp, "TimeLimit", 60*60)
	set_optimizer_attribute(smp, "OutputFlag", 0)

	#Variables
	x = Dict()
	for i in orders, a in orderArcSet[i]
	    global x[i,a] = @variable(smp, lower_bound = 0) #, upper_bound = 1) <-- UB implied by order origin/dest con & messes up reduced cost sign
	    set_name(x[i,a], string("x[",i,",",a,"]")) 
	end
	@variable(smp, y[A_hasdriver] >= 0)
	@variable(smp, z[l = 1:numlocs, s = 1:numshifts, f = 1:numfragments[l,s]] >= 0)
	#@variable(smp, 0 <= z[d = drivers, homeArcSet[d]] <= 1)
	@variable(smp, w[a in A_space] >= 0)
	@variable(smp, ordtime[orders])

	#Objective
	@objective(smp, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderArcSet[i]) for i in orders) + sum(c[a]*y[a] for a in A_hasdriver) + sum(u[a]*w[a] for a in A_space) )
	
	#Order constraints
	@constraint(smp, orderFlowBalance[i = orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))], sum(x[i,a] for a in A_minus_i[i,n]) - sum(x[i,a] for a in A_plus_i[i,n]) == 0)
	@constraint(smp, arriveDestin[i = orders], sum(sum(x[i,a] for a in A_minus_i[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in Destination[i]) == 1)
	@constraint(smp, departOrigin[i = orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), A_plus_i[i,n])) for n in Origin[i]) == 1)
	for i in setdiff(orders, ordersinprogress)
		extendedorderarc = extendedarcs[last(Origin[i]), last(Destination[i])]
		if extendedorderarc in orderArcSet[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end
	
	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(orders, ordersinprogress), n in Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
		if a in orderArcSet[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(smp, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in A_minus_i[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(smp, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_plus_hd[n]) == m_0[n])
	@constraint(smp, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) >= m_end[n])
	@constraint(smp, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n])- sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n],dummyarc)) for i in orders) - sum(y[a] for a in A_plus_hd[n]) == 0)

	#Linking constraints
	@constraint(smp, driverAvailability[a in A_space], sum(sum(sum(z[l,s,f] for f in fragmentscontaining[l,s,a]) for s in 1:numshifts) for l in 1:numlocs) == w[a] )
	for i in orders, a in orderArcSet_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in A_hasdriver_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(smp, driverStartingLocs[l in 1:numlocs, s in 1:numshifts], sum(sum(z[l,s,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == length(driversets[l,s]))
	@constraint(smp, driverFlowBalance[l in 1:numlocs, s in 1:numshifts, n in N_flow_ls[l,s]], sum(z[l,s,f] for f in F_minus_ls[l,s,n]) - sum(z[l,s,f] for f in F_plus_ls[l,s,n]) == 0)

	#Return named tuple of constraints needed for column generation
	smpconstraints = (con_orderFlowBalance = orderFlowBalance,
		con_departOrigin = departOrigin,
		con_arriveDestin = arriveDestin,
		con_initialTrucks = initialTrucks,
		con_finalTrucks = finalTrucks,
		con_truckFlowBalance = truckFlowBalance,
		con_driverAvailability = driverAvailability,
		con_deliveryTime = deliveryTime
		)

	return smp, x, y, z, w, smpconstraints

end

#----------------------------------------------------------------------------------------#

function getarcsetforsubproblem(i, arcredcosts, orderArcSet)

	#Create copies of all important sets (so they can be modified for the shortest path problem)
	arcsSP, arcLookupSP, nodesLookupSP = deepcopy(arcs), deepcopy(arcLookup), deepcopy(nodesLookup)

	arclistSP = deepcopy(setdiff(orderArcSet[i], dummyarc))
	cSP = Dict()
	for a in arclistSP
		cSP[a] = arcredcosts[i,a]
	end
	A_minusSP = Dict()
	for n in 1:extendednumnodes
		A_minusSP[n] = copy(intersect(A_minus_i[i,n], setdiff(orderArcSet[i], dummyarc))) 
		if dummyarc in A_minusSP[n]
			println("Here's the culprit!! --> $n")
		end
	end

	#Create dummy origin and destination nodes + arcs 
	dummyorig = extendednumnodes + 1
	dummydest = extendednumnodes + 2
	numnodesSP = extendednumnodes + 2
	nodesLookupSP[dummyorig], nodesLookupSP[dummydest] = (-1,-1), (-1,99999)
	A_minusSP[dummyorig] = []
	A_minusSP[dummydest] = []

	#Add arcs from dummy origin to all origin nodes in pickup window
	index = extendednumarcs + 2
	for j in 1:length(Origin[i])
		n = Origin[i][j]
		arcsSP[dummyorig, n] = index
		arcLookupSP[index] = (dummyorig, n)
		cSP[index] = 0
		pushfirst!(arclistSP, index)
		push!(A_minusSP[n], index)
		index += 1
	end

	#Add arcs from all destination nodes in delivery window to the dummy destination
	for n in Destination[i]
		arcsSP[n, dummydest] = index
		arcLookupSP[index] = (n, dummydest) 
		cSP[index] = 0
		push!(arclistSP, index)
		push!(A_minusSP[dummydest], index)
		index += 1
	end
	
	#Remove arcs that are not feasible because they are longer than the driver shift
	for a in 1:numarcs
		arctime = nodesLookup[arcLookup[a][2]][2] - nodesLookup[arcLookup[a][1]][2]
		if arctime > shiftlength
			remove!(arclistSP, a)
			remove!(A_minusSP[arcLookup[a][2]], a)
		end
	end

	for n in 1:extendednumnodes
		if dummyarc in A_minusSP[n]
			println("Later... Here's the culprit!! --> $i, $n")
		end
	end

	return numnodesSP, arclistSP, nodesLookupSP, arcLookupSP, cSP, dummyorig, dummydest, A_minusSP

end

#----------------------------------------------------------------------------------------#

function mvgsubproblem_preprocess(orderArcSet)

	#Create copies of all important sets (so they can be modified for the shortest path problem)
	nodesLookupSP = deepcopy(nodesLookup)
	arcsSP, arcLookupSP, nodesLookupSP = Dict(), Dict(), Dict()
	for i in orders
		arcsSP[i], arcLookupSP[i] = deepcopy(arcs), deepcopy(arcLookup)
	end

	#Create arclists
	arclistSP = Dict()
	for i in orders
		arclistSP[i] = deepcopy(setdiff(orderArcSet[i], dummyarc))
	end

	#Create dummy origin and destination nodes + arcs 
	dummyorig = extendednumnodes + 1
	dummydest = extendednumnodes + 2
	numnodesSP = extendednumnodes + 2
	nodesLookupSP[dummyorig], nodesLookupSP[dummydest] = (-1,-1), (-1,99999)

	#Add dummy arcs
	for i in orders
		index = extendednumarcs + 2

		#Add arcs from dummy origin to all origin nodes in pickup window
		for j in 1:length(Origin[i])
			n = Origin[i][j]
			arcsSP[i][dummyorig, n] = index
			arcLookupSP[i][index] = (dummyorig, n)
			pushfirst!(arclistSP[i], index)
			index += 1
		end

		#Add arcs from all destination nodes in delivery window to the dummy destination
		for n in Destination[i]
			arcsSP[i][n, dummydest] = index
			arcLookupSP[i][index] = (n, dummydest) 
			push!(arclistSP[i], index)
			index += 1
		end
	end

	#Remove arcs that are not feasible because they are longer than the driver shift
	for i in orders, a in 1:numarcs
		arctime = nodesLookup[arcLookup[a][2]][2] - nodesLookup[arcLookup[a][1]][2]
		if arctime > shiftlength
			remove!(arclistSP[i], a)
		end
	end

	return numnodesSP, arclistSP, nodesLookupSP, arcLookupSP, dummyorig, dummydest

end

#----------------------------------------------------------------------------------------#

function mvgsubproblem_fast(i, arcredcosts, orderArcSet, numnodesSP, arclistSP, nodesLookupSP, arcLookupSP, dummyorig, dummydest)

	cSP = Dict()
	for a in arclistSP
		try
			cSP[a] = arcredcosts[i,a]
		catch
			cSP[a] = 0 #dummy orig and dest
		end	
	end

	#Initialize shortest path algorithm 
	currdistance = repeat([999999999.0],outer=[numnodesSP])
	currdistance[dummyorig] = 0
	prevnode, prevarc = zeros(numnodesSP), zeros(numnodesSP)
	orderstarttime = nodesLookup[Origin[i][1]][2]
	maxspiter = (horizon - orderstarttime)/tstep + 3
	
	#Loop over time-space arcs in order of start time
	for a in arclistSP 
		n_end, n_start = arcLookupSP[a][2], arcLookupSP[a][1]
		if currdistance[n_end] > currdistance[n_start] + cSP[a] + 0.000001
			currdistance[n_end] = currdistance[n_start] + cSP[a]
			prevnode[n_end] = n_start
			prevarc[n_end] = a
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

	return currdistance[dummydest], shortestpathnodes, shortestpatharcs

end

#----------------------------------------------------------------------------------------#

function mvgsubproblem(i, arcredcosts, orderArcSet)

	spstart = time()
	numnodesSP, arclistSP, nodesLookupSP, arcLookupSP, cSP, dummyorig, dummydest = getarcsetforsubproblem(i, arcredcosts, orderArcSet)
	time1 = spstart - time()

	#Initialize shortest path algorithm 
	currdistance = repeat([999999999.0],outer=[numnodesSP])
	currdistance[dummyorig] = 0
	prevnode, prevarc = zeros(numnodesSP), zeros(numnodesSP)
	orderstarttime = nodesLookup[Origin[i][1]][2]
	maxspiter = (horizon - orderstarttime)/tstep + 3
	time2 = time1 - time()
	
	#Bellman Ford loop
	#for iteration in 1:maxspiter #Max path length
	for a in arclistSP #1:numarcsSP
		#println("arcLookupSP[$a] = ", arcLookupSP[a])
		n_end, n_start = arcLookupSP[a][2], arcLookupSP[a][1]
		if currdistance[n_end] > currdistance[n_start] + cSP[a] + 0.000001
			currdistance[n_end] = currdistance[n_start] + cSP[a]
			prevnode[n_end] = n_start
			prevarc[n_end] = a
		end
	end
	#end
	time3 = time2 - time()

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
	time4 = time3 - time()

	return currdistance[dummydest], shortestpathnodes, shortestpatharcs, time1, time2, time3, time4

end

#----------------------------------------------------------------------------------------#

function mvgsubproblem_bernardo(i, arcredcosts, orderArcSet)

	#Fill in later to see if we can speed things up
	numnodesSP, arclistSP, nodesLookupSP, arcLookupSP, cSP, dummyorig, dummydest, A_minusSP = getarcsetforsubproblem(i, arcredcosts, orderArcSet)

	#Initialize shortest path algorithm 
	currdistance = repeat([999999999.0],outer=[numnodesSP])
	currdistance[dummyorig] = 0
	prevnode, prevarc = zeros(numnodesSP), zeros(numnodesSP)
	orderstarttime = nodesLookup[Origin[i][1]][2]
	maxspiter = (horizon - orderstarttime)/tstep + 3

	#sort!(arclistSP, by = x -> nodesLookupSP[arcLookupSP[x][1]][2])
	for a in arclistSP #1:numarcsSP
		#println("arcLookupSP[$a] = ", arcLookupSP[a])
		n_end, n_start = arcLookupSP[a][2], arcLookupSP[a][1]
		if currdistance[n_end] > currdistance[n_start] + cSP[a] + 0.000001
			currdistance[n_end] = currdistance[n_start] + cSP[a]
			prevnode[n_end] = n_start
			prevarc[n_end] = a
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

	return currdistance[dummydest], shortestpathnodes, shortestpatharcs	

end

#----------------------------------------------------------------------------------------#

function findarcvariablereducedcosts(z, alpha, beta, gamma, theta, nu, mu, xi, psi, bpRestrictedArcSet)

	arcredcosts = Dict()
	for i in orders, a in 1:extendednumarcs
		arcredcosts[i,a] = c[a] 
	end

	#Calculate reduced costs for each arc
	for i in orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i])), a in A_minus[n]
		arcredcosts[i,a] -= alpha[i,n]
	end
	for i in orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i])), a in A_plus[n]
		arcredcosts[i,a] += alpha[i,n] 
	end

	for i in orders, n in Origin[i], a in intersect(union(A_space, dummyarc), A_plus[n])
		arcredcosts[i,a] -= beta[i] 
	end
	for i in setdiff(orders, ordersinprogress)
		a = extendedarcs[last(Origin[i]), last(Destination[i])]
		arcredcosts[i,extendedarcs[last(Origin[i]), last(Destination[i])]] -= beta[i] 
	end
	
	for i in intersect(orders, ordersinprogress), n in Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
		arcredcosts[i,a] -= beta[i] 
	end
	
	for i in orders, n in Destination[i], a in A_minus[n]
		if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))
			arcredcosts[i,a] -= gamma[i]
		end
	end
	
	for n in N_0, i in orders, a in setdiff(A_plus[n], dummyarc)  
		arcredcosts[i,a] -= theta[n]
	end
	for n in N_end, a in A_minus[n], i in orders
		arcredcosts[i,a] -= nu[n]
	end
	for n in N_flow_t, a in A_minus[n], i in orders
		arcredcosts[i,a] -= mu[n]
	end
	for n in N_flow_t, a in A_plus[n], i in orders
		arcredcosts[i,a] += mu[n]
	end
	
	for i in orders, n in Destination[i], a in A_minus[n]
		arcredcosts[i,a] += arcfinishtime[a] * psi[i]
	end
	for i in orders, a in A_space
		arcredcosts[i,a] += xi[a]
	end

	return arcredcosts

end

#----------------------------------------------------------------------------------------#

function formatduals(alpha_stg, beta_stg, gamma_stg, theta_stg, nu_stg, mu_stg, xi_stg, psi_stg)

	#Alpha
	alpha = Dict()
	alpha_row, alpha_col, alpha_val = Dict(),Dict(),Dict()
	for i in orders
		alpha_row[i] = []
		alpha_col[i] = []
		alpha_val[i] = []
	end
	for (i,n) in eachindex(alpha_stg)
		push!(alpha_row[i], n)
		push!(alpha_col[i], 1)
		push!(alpha_val[i], alpha_stg[i,n])
	end
	for i in orders
		alpha[i] = sparse(alpha_row[i], alpha_col[i], alpha_val[i])
	end

	#Beta
	beta_row, beta_col, beta_val = [], [], []
	for i in orders
		push!(beta_row, i)
		push!(beta_col, 1)
		push!(beta_val, beta_stg[i])
	end
	beta = sparse(beta_row, beta_col, beta_val)

	#Gamma
	gamma_row, gamma_col, gamma_val = [], [], []
	for i in orders
		push!(gamma_row, i)
		push!(gamma_col, 1)
		push!(gamma_val, gamma_stg[i])
	end
	gamma = sparse(gamma_row, gamma_col, gamma_val)

	#Theta
	theta_row, theta_col, theta_val = [n for n in N_0], [1 for n in N_0], [theta_stg[n] for n in N_0]
	theta = sparse(theta_row, theta_col, theta_val)

	#Nu
	nu_row, nu_col, nu_val = [n for n in N_end], [1 for n in N_end], [nu_stg[n] for n in N_end]
	nu = sparse(nu_row, nu_col, nu_val)

	#Mu
	mu_row, mu_col, mu_val = [n for n in N_flow_t], [1 for n in N_flow_t], [mu_stg[n] for n in N_flow_t]
	mu = sparse(mu_row, mu_col, mu_val)

	#Xi
	xi_row, xi_col, xi_val = [a for a in A_space], [1 for a in A_space], [xi_stg[a] for a in A_space]
	xi = sparse(xi_row, xi_col, xi_val)

	#Psi
	psi_row, psi_col, psi_val = [i for i in orders], [1 for i in orders], [psi_stg[i] for i in orders]
	psi = sparse(psi_row, psi_col, psi_val)

	return alpha, beta, gamma, theta, nu, mu, xi, psi
		
end

#----------------------------------------------------------------------------------------#

function findarcvariablereducedcosts_fast(z, alpha, beta, gamma, theta, nu, mu, xi, psi, bpRestrictedArcSet)

	#Initialize
	arcredcosts = Dict()
	for i in orders, a in 1:extendednumarcs
		arcredcosts[i,a] = c[a]
	end

	#Add constraint dual contributions
	#Alpha
	for i in orders
		product_alpha = xcoeffs_alpha[i] * alpha[i] 
		for a in rowvals(product_alpha)
			arcredcosts[i,a] -= product_alpha[a,1]
		end
	end

	#Beta
	product_beta = dropzeros(xcoeffs_beta * beta)
	nza, nzi, nzv = findnz(product_beta)
	for index in 1:length(nza)
		arcredcosts[nzi[index], nza[index]] -= nzv[index]
	end

	#Gamma
	product_gamma = dropzeros(xcoeffs_gamma * gamma)
	nza, nzi, nzv = findnz(product_gamma)
	for index in 1:length(nza)
		arcredcosts[nzi[index], nza[index]] -= nzv[index]
	end

	#Theta
	for i in orders
		nza, nzn, nzv = findnz(xcoeffs_theta[i])
		for index in 1:length(nza)
			arcredcosts[i, nza[index]] -= nzv[index] * theta[nzn[index]]
		end
	end

	#Nu
	for i in orders
		nza, nzn, nzv = findnz(xcoeffs_nu[i])
		for index in 1:length(nza)
			arcredcosts[i, nza[index]] -= nzv[index] * nu[nzn[index]]
		end
	end

	#Mu
	for i in orders
		nza, nzn, nzv = findnz(xcoeffs_mu[i])
		for index in 1:length(nza)
			arcredcosts[i, nza[index]] -= nzv[index] * mu[nzn[index]]
		end
	end

	#Xi
	nza, nzi, nzv = findnz(xcoeffs_xi)
	for index in 1:length(nza)
		arcredcosts[nzi[index], nza[index]] -= nzv[index] * xi[nza[index], 1]
	end

	#Psi
	nza, nzi, nzv = findnz(xcoeffs_psi)
	for index in 1:length(nza)
		arcredcosts[nzi[index], nza[index]] -= nzv[index] * psi[nzi[index]]
	end

	#Incorporate driver info into the reduced cost consideration
	if newreducedcost_flag == 1
		driverredcost = reduced_cost.(z)
		bestdrivercost = []
		for a in 1:numarcs
			if (length(availableDrivers[a]) > 0) & (a in A_space)
				push!(bestdrivercost, max(0,minimum(driverredcost[d,a] for d in availableDrivers[a])))
			else
				push!(bestdrivercost, 0)
			end
		end

		for i in orders, a in 1:numarcs
			arcredcosts[i,a] += bestdrivercost[a]
		end
	end

	return arcredcosts

end

#----------------------------------------------------------------------------------------#

function strengthenreducedcosts(arcredcosts, z, orderArcSet_full)

	driverredcost = reduced_cost.(z)

	bestdriver_rc = Dict()
	for a in 1:numarcs
		listoffragments = []
		for l in 1:numlocs, s in 1:numshifts, f in fragmentscontaining[l,s,a]
			push!(listoffragments, (l,s,f))
		end
		if listoffragments != []
			bestdriver_rc[a] = max(0, minimum(driverredcost[frag] for frag in listoffragments))
		else
			bestdriver_rc[a] = 0
		end
	end

	for i in orders, a in 1:numarcs
		arcredcosts[i,a] += bestdriver_rc[a]
	end

	return arcredcosts

end

#----------------------------------------------------------------------------------------#

function fragmentmultivariablegeneration!(orderArcSet, orderArcSet_space, A_plus_i, A_minus_i, bpRestrictedArcSet, homeArcSet, A_minus_d, A_plus_d, availableDrivers, A_minus_i_full, A_plus_i_full)

	#Build sparse master problem
	model, x, y, z, w, smpconstraints = sparsemasterproblem(orderArcSet, orderArcSet_space, A_plus_i, A_minus_i, homeArcSet, A_minus_d, A_plus_d, availableDrivers)

	#------------------------------------------------------#

	#Initialize column generation 
	cg_iter = 1
	smpobjectives, smptimes, pptimes, pptimes_par, lowerbounds = [], [], [], [], []
	listlength = convert(Int64, ceil(length(orders)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)
	currminredcost = -1000000000

	#------------------------------------------------------#

	#Pre-processing - constraint matrix + subproblem
	#xcoeffs_alpha, xcoeffs_beta, xcoeffs_gamma, xcoeffs_theta, xcoeffs_nu, xcoeffs_mu, xcoeffs_xi, xcoeffs_psi = findxcoefficientmatrix(A_minus_i_full, A_plus_i_full) 
	numnodesSP, arclistSP, nodesLookupSP, arcLookupSP, dummyorig, dummydest = mvgsubproblem_preprocess(bpRestrictedArcSet)

	#------------------------------------------------------#

	while cg_iter <= 100000

		#-----------SOLVE SMP-----------#

		#println("-------- ITERATION $cg_iter --------")

		status = optimize!(model)
		if termination_status(model) != MOI.OPTIMAL
			println(termination_status(model))
			return 100000000, model, x, y, z, w, orderArcSet
		end
		smpobj, smptime = objective_value(model), solve_time(model)
		
		push!(smpobjectives, copy(smpobj))
		push!(smptimes, copy(smptime))	
		if printstatements == 1
			println("Solved SMP with objective = ", smpobj, " in iteration $cg_iter (", sum(length(orderArcSet[i]) for i in orders), " arcs)")
			#println("SMP time = ", smptime)
		end

		#------------SNAPSHOT-----------#
		starttime1 = time()
		orderArcSet_snapshot = deepcopy(orderArcSet)
		fulltime1 = time() - starttime1
		#println("Snapshot time = ", fulltime1)

		#---------RED COST CALC---------#
	
		#Get constraint dual values
		alpha = dual.(smpconstraints.con_orderFlowBalance)
		beta = dual.(smpconstraints.con_departOrigin)
		gamma = dual.(smpconstraints.con_arriveDestin)
		theta = dual.(smpconstraints.con_initialTrucks)
		nu = dual.(smpconstraints.con_finalTrucks)
		mu = dual.(smpconstraints.con_truckFlowBalance)
		xi = dual.(smpconstraints.con_driverAvailability)
		psi = dual.(smpconstraints.con_deliveryTime)

		#Get reduced costs
		starttime3 = time()
		arcredcosts = findarcvariablereducedcosts(z, alpha, beta, gamma, theta, nu, mu, xi, psi, bpRestrictedArcSet)
		if newreducedcost_flag == 1
			arcredcosts = strengthenreducedcosts(arcredcosts, z, orderArcSet_full)
		end
		fulltime3 = time() - starttime3
		#println("Arc reduced costs time = ", fulltime3)

		#starttime3_new = time()
		#alpha_new, beta_new, gamma_new, theta_new, nu_new, mu_new, xi_new, psi_new = formatduals(alpha, beta, gamma, theta, nu, mu, xi, psi)
		#findarcvariablereducedcosts_fast(z, alpha_new, beta_new, gamma_new, theta_new, nu_new, mu_new, xi_new, psi_new, bpRestrictedArcSet)
		#fulltime3_new = time() - starttime3_new
		#println("Arc reduced costs time (fast) = ", fulltime3_new)

		#------------SUBPROBLEM------------#

		#Run shortest path for each order to find new arcs
		starttime4 = time()
		
		dptimelist, dptimelist_ber = [], []
		addarcs, minreducedcosts = [], []
		for i in orders
			starttime = time()
			if onearcatatime_flag == 0	
				minreducedcost, shortestpathnodes, shortestpatharcs = mvgsubproblem_fast(i, arcredcosts, bpRestrictedArcSet, numnodesSP, arclistSP[i], nodesLookupSP, arcLookupSP[i], dummyorig, dummydest)
			elseif onearcatatime_flag == 1
				mostnegativearc = setdiff(orderArcSet_full[i], dummyarc)[argmin([arcredcosts[i,a] for a in setdiff(orderArcSet_full[i], dummyarc)])]
				minreducedcost = arcredcosts[i, mostnegativearc]
				shortestpatharcs = [mostnegativearc]
				shortestpathnodes = []
			end
			endtime = time() - starttime

			push!(dptimelist, endtime)
			push!(minreducedcosts, minreducedcost)

			if minreducedcost < -0.0001
				for a in shortestpatharcs
					push!(addarcs, (i,a))
				end				
			end
		end
		push!(pptimes, sum(dptimelist))
		#push!(pptimes_ber, sum(dptimelist_ber))
		#println("Us = ", sum(dptimelist), " vs. bernardo = ", sum(dptimelist_ber))

		#"Parallelize"
		shuffleddptimes = shuffle_partition(length(orders))
		dptimelistsums = [sum(dptimelist[j_inner] for j_inner in shuffleddptimes[j_outer]) for j_outer in 1:length(shuffleddptimes)]
		push!(pptimes_par, maximum(dptimelistsums))
		try
			push!(lowerbounds, smpobj + sum(minreducedcosts[k] for k in 1:length(minreducedcosts) if minreducedcosts[k] <= -0.000001))
		catch
			push!(lowerbounds, smpobj)
		end
		#println("PP time (parallel)= ", maximum(dptimelistsums))

		fulltime4 = time() - starttime4
		#println("Full SP time = ", fulltime4)

		#------COUNT ARCS AND PATHS-----#

		if saveconvergencedata_flag >= 0
			totalorderarcs = sum(length(orderArcSet[i]) for i in orders)
			totalorderpaths = sum([findallpaths(A_plus_i, i) for i in orders])
			maximprove = minimum(minreducedcosts) * sum(sum(value(x[i,a]) for a in orderArcSet[i]) for i in orders)
			write_cg_conv(convergencedatafilename, cg_iter, maximprove, totalorderarcs, totalorderpaths, smpobj)
		end
	
		#-------ADD NEW VARIABLES-------#

		starttime5 = time()
		#Add new arcs to order arc sets
		for (i,a) in addarcs
			if !(a in orderArcSet[i])
				push!(orderArcSet[i], a)
				if a in A_space
					push!(orderArcSet_space[i], a)
				end
				n_plus = arcLookup[a][1]
				n_minus = arcLookup[a][2]
				push!(A_plus_i[i,n_plus], a)
				push!(A_minus_i[i,n_minus], a)
			end
		end

		fulltime5 = time() - starttime5
		#println("Add new arcs time = ", fulltime5)

		starttime6 = time()

		#Add new arcs to model as x-variables
		for i in orders, a in setdiff(orderArcSet[i], orderArcSet_snapshot[i])

			#Create a new variable for the arc
			global x[i,a] = @variable(model, lower_bound = 0) #, upper_bound = 1)
			set_name(x[i,a], string("x[",i,",",a,"]")) 

			#Add to the objective
			set_objective_function(model, objective_function(model) + c[a]*x[i,a])

			#Find the entering and exiting nodes for arc a (ex. n_plus is the node for which a belongs to A_plus[n])
			n_plus = arcLookup[a][1]
			n_minus = arcLookup[a][2]

			#Add new variable to order constraints
			if n_minus in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))
				set_normalized_coefficient(smpconstraints.con_orderFlowBalance[i,n_minus], x[i,a], 1.0)
			end
			if n_plus in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))
				set_normalized_coefficient(smpconstraints.con_orderFlowBalance[i,n_plus], x[i,a], -1.0)
			end
			if (n_plus in Origin[i]) & (a in A_space)
				set_normalized_coefficient(smpconstraints.con_departOrigin[i], x[i,a], 1.0)
			end
			if (n_plus in Origin[i]) & (i in intersect(orders, ordersinprogress)) & !(a in union(A_space, dummyarc))
				set_normalized_coefficient(smpconstraints.con_departOrigin[i], x[i,a], 1.0)
			end
			if (inprogressdummyarc_flag == 0) && !(i in ordersinprogress) && (a == extendedarcs[last(Origin[i]), last(Destination[i])]) 
				set_normalized_coefficient(smpconstraints.con_departOrigin[i], x[i,a], 1.0)
			end
			if (n_minus in Destination[i]) & ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))
				set_normalized_coefficient(smpconstraints.con_arriveDestin[i], x[i,a], 1.0)
			end

			#Add new variable to order delivery constraints
			if (n_minus in Destination[i]) & (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1])
				set_normalized_coefficient(smpconstraints.con_deliveryTime[i], x[i,a], -arcfinishtime[a] )
			end
		
			#Add new variable to initial and final trucks constraints
			if n_plus in N_0
				set_normalized_coefficient(smpconstraints.con_initialTrucks[n_plus], x[i,a], 1.0)
			end
			if n_minus in N_end
				set_normalized_coefficient(smpconstraints.con_finalTrucks[n_minus], x[i,a], 1.0)
			end

			#Add new variable to truck flow balance constraints
			if n_minus in N_flow_t
				set_normalized_coefficient(smpconstraints.con_truckFlowBalance[n_minus], x[i,a], 1.0 )
			end
			if n_plus in N_flow_t
				set_normalized_coefficient(smpconstraints.con_truckFlowBalance[n_plus], x[i,a], -1.0 )
			end

			#Add new variable to driver availability constraints
			if a in A_space
				set_normalized_coefficient(smpconstraints.con_driverAvailability[a], x[i,a], -1.0)
			end

		end

		fulltime6 = time() - starttime6
		#println("Add vars to SMP time = ", fulltime6)

		#----------TERMINATION----------#

		starttime7 = time()
		if minimum(minreducedcosts) >= -0.0001
			if printstatements == 1
				println("NO NEGATIVE REDUCED COSTS FOUND!")	
			end

			break

		end

		#------------ITERATE------------#

		cg_iter += 1

		fulltime7 = time() - starttime7
		#println("Transition time = ", fulltime7)

	end

	optimize!(model)
	smpobj = objective_value(model)

	#-------------------------------#	

	df = DataFrame(experiment_id = [experiment_id for k in 1:length(smptimes)],
			instance = [ex for j in 1:length(smptimes)],
			lambda = [lambda for j in 1:length(smptimes)],
			horizon = [horizon for j in 1:length(smptimes)],
			tstep = [tstep for j in 1:length(smptimes)],
			week = [weekstart for j in 1:length(smptimes)],
			driverfactor = [driverfactor for j in 1:length(smptimes)],
			numdrivers = [length(drivers) for j in 1:length(smptimes)],
			method = [solutionmethod for j in 1:length(smptimes)],
			k = [k for j in 1:length(smptimes)],
			iteration = [j for j in 1:length(smptimes)],
            objective = smpobjectives,
            milesobj = [0 for j in 1:length(smptimes)],
            delayobj = [0 for j in 1:length(smptimes)],
            bound = lowerbounds,
            rmptime = smptimes,
            pptime = pptimes,
            pptime_par = pptimes_par,
            iptime = [0 for j in 1:length(smptimes)],
			mvgtime_full = [0 for j in 1:length(smptimes)]
           )

	CSV.write(resultsfilename, df)

	#-------------------------------#

	return smpobj, model, x, y, z, w, orderArcSet, sum(smptimes), sum(pptimes), sum(pptimes_par)

end

#----------------------------------------------------------------------------------------#

function applyarcrestrictions(orderArcSet_rest, orderArcSet_space_full, A_plus_i_full, A_minus_i_full)

	orderArcSet_space_rest, A_plus_i_rest, A_minus_i_rest = Dict(), Dict(), Dict()

	for i in orders
		orderArcSet_space_rest[i] = intersect(orderArcSet_rest[i], orderArcSet_space_full[i])
		for n in 1:extendednumnodes
			A_plus_i_rest[i,n] = intersect(orderArcSet_rest[i], A_plus_i_full[i,n])
			A_minus_i_rest[i,n] = intersect(orderArcSet_rest[i], A_minus_i_full[i,n])
		end
	end	

	return orderArcSet_space_rest, A_plus_i_rest, A_minus_i_rest

end