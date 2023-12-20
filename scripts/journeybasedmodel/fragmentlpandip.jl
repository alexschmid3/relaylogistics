
function createdriversets()
		
	#Find the effective shifts 
	#i.e. For a three day horizon starting on a Monday, a MTWTF driver looks the same as a MTWR+Su driver
	effectiveshifts = unique!([T_off[s] for s in 1:numshifts])
	numeffshifts = length(effectiveshifts)
	effshift, shiftsincluded = [], [[] for es in 1:numeffshifts]
	for s in 1:numshifts
		es = findfirst(x->x==T_off[s], effectiveshifts)
		push!(effshift, es)
		push!(shiftsincluded[es], s)
	end

	#Create driver sets
	driversets = Dict()
	for l in 1:numlocs, s in 1:numeffshifts
		driversets[l,s] = []
	end
	for d in drivers
		push!(driversets[driverHomeLocs[d], effshift[drivershift[d]]], d)
	end

	driverSetStartNodes = Dict()
	for l in 1:numlocs, s in 1:numeffshifts
		driverSetStartNodes[l,s] = nodes[l,0]
	end

	return driversets, driverSetStartNodes, numeffshifts, effshift, shiftsincluded

end

#----------------------------------------------------------------------------------------#

function addnextsegment(currentfrags, completedfrags, adjacentlocs, deadline, homeloc)
	
	currfragnum = length(currentfrags)

	for f in 1:currfragnum
		frag = popfirst!(currentfrags)
		currloc, currtime = last(frag)
		if mod(currtime, 24) < shiftlength #Working hours
			for a in union(adjacentlocs[currloc], currloc)
				if (a == currloc) & (currtime + tstep <= deadline - traveltimebetweenlocs_rdd[a, homeloc]) 
					newfrag = deepcopy(frag)
					push!(newfrag, (a, currtime + tstep))
					push!(currentfrags, newfrag)
					if currtime + tstep == deadline
						push!(completedfrags, newfrag)
					end
				elseif (a != currloc) && (mod(currtime + arcLength[currloc, a],24) <= shiftlength) && (currtime + arcLength[currloc, a] <= deadline - traveltimebetweenlocs_rdd[a, homeloc]) 
					newfrag = deepcopy(frag)
					push!(newfrag, (a, currtime + arcLength[currloc, a]))
					push!(currentfrags, newfrag)
					if currtime + arcLength[currloc, a] == deadline
						push!(completedfrags, newfrag)
					end
				end
			end
		elseif currloc != homeloc
			newfrag = deepcopy(frag)
			push!(newfrag, (currloc, currtime + tstep))
			push!(currentfrags, newfrag)
			if currtime + tstep == deadline
				push!(completedfrags, newfrag)
			end
		end
	end
	
	currentfrags = setdiff(currentfrags, completedfrags)

	return currentfrags, completedfrags

end

#----------------------------------------------------------------------------------------#

function findbasefragments(homeloc, maxnightsaway, adjacentlocs)
		
	basefragments = []

	for deadline in shiftlength:24:shiftlength+24*maxnightsaway
		currentfrags = [[(homeloc, 0)]] 
		completedfrags = []
		
		#Loop to enumerate paths
		while currentfrags != []
			currentfrags, completedfrags = addnextsegment(currentfrags, completedfrags, adjacentlocs, deadline, homeloc)
	 	end

	 	#Add the paths of length "deadline" to the list of fragments 
	 	basefragments = union(basefragments, completedfrags)
	end

	return basefragments

end

#----------------------------------------------------------------------------------------#

function createfragmentsets(maxnightsaway)

	driversets, driverSetStartNodes, numeffshifts, effshift, shiftsincluded = createdriversets()

	adjacentlocs = [[] for l in 1:numlocs]
	for pa in prearcs
		if pa[3] <= shiftlength
			push!(adjacentlocs[pa[1]], pa[2])
		end
	end

	numfragments = Dict()
	fragmentscontaining = Dict()
	F_plus_ls, F_minus_ls = Dict(), Dict()
	for l in 1:numlocs, s in 1:numeffshifts
		numfragments[l,s] = 0
	end
	for l in 1:numlocs, s in 1:numeffshifts, a in 1:numarcs
		fragmentscontaining[l,s,a] = []
	end
	for l in 1:numlocs, s in 1:numeffshifts, n in 1:numnodes
		F_plus_ls[l,s,n] = []
		F_minus_ls[l,s,n] = []
	end

	for l in 1:numlocs
		basefragments = findbasefragments(l, maxnightsaway, adjacentlocs)
		for s in 1:numeffshifts
			#println("------ $l, $s ------")
			fragindex = 1
			if driversets[l,s] != []
				d_ex = driversets[l,s][1]
				#Dummy fragment to stay home during off hours at time 0
				if !(0 in T_on_0[d_ex])
					frag = [(l, t) for t in 0:tstep:T_on_0[d_ex][1]]

					#Create official timed fragment 
					currnode = nodes[frag[1][1], frag[1][2]]
					push!(F_plus_ls[l,s,currnode], fragindex)

					for n in 2:length(frag)
						newnode = nodes[frag[n][1], frag[n][2]]
						newarc = arcs[currnode, newnode]
						push!(fragmentscontaining[l,s,newarc], fragindex)
						currnode = newnode
					end

					push!(F_minus_ls[l,s,currnode], fragindex)

					#println(fragindex, " ==> ", frag, " @ t = 0 (pre-shift)")

					fragindex += 1
					numfragments[l,s] += 1

				end	

				for t in [t2 for t2 in T_on_0[d_ex] if t2 <= horizon] #All possible start times of plans
					for frag in basefragments
						fragmentendtime = t + last(frag)[2]
						
						if (fragmentendtime - shiftlength in T_on_0[d_ex]) || (fragmentendtime - shiftlength in T_on_0[d_ex])
							#Create official timed fragment 
							currnode = nodes[frag[1][1], frag[1][2] + t]
							push!(F_plus_ls[l,s,currnode], fragindex)

							for n in 2:length(frag)
								if frag[n][2] + t <= horizon
									newnode = nodes[frag[n][1], frag[n][2] + t]
									newarc = arcs[currnode, newnode]
									push!(fragmentscontaining[l,s,newarc], fragindex)
									currnode = newnode
								end
							end

							#Add the off hours rest period to the end of each fragment to connect to the start of the next fragment
							for trest in fragmentendtime+tstep:tstep:min(horizon, fragmentendtime+24-shiftlength)
								newnode = nodes[l, trest]
								newarc = arcs[currnode, newnode]
								push!(fragmentscontaining[l,s,newarc], fragindex)
								currnode = newnode
							end

							push!(F_minus_ls[l,s,currnode], fragindex)

							#println(fragindex, " ==> ", frag, " @ t = $t")

							fragindex += 1
							numfragments[l,s] += 1
						end
					end
				end
			end
		end
	end

	N_flow_ls = Dict()
	for l in 1:numlocs, s in 1:numeffshifts
		N_flow_ls[l,s] = []
		for n in setdiff(1:numnodes, union(driverSetStartNodes[l,s], N_end))
			#if union(F_plus_ls[l,s,n], F_minus_ls[l,s,n]) != []
			push!(N_flow_ls[l,s], n)
			#end
		end
	end

	return driversets, driverSetStartNodes, numfragments, fragmentscontaining, F_plus_ls, F_minus_ls, N_flow_ls, numeffshifts, effshift, shiftsincluded

end 

#----------------------------------------------------------------------------------------#

function solvefragmentlp(opt_gap, orderArcSet, orderArcSet_space, A_plus_i, A_minus_i, numeffshifts)

	lp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(lp, "TimeLimit", 60*60*40)
	set_optimizer_attribute(lp, "OutputFlag", 0)
	set_optimizer_attribute(lp, "MIPGap", opt_gap)

	#Variables
	@variable(lp, 0 <= x[i in orders, a in orderArcSet[i]] <= 1)
	@variable(lp, y[A_hasdriver] >= 0)
	@variable(lp, z[l = 1:numlocs, s = 1:numeffshifts, f = 1:numfragments[l,s]] >= 0)
	@variable(lp, w[a in A_space] >= 0)
	@variable(lp, ordtime[orders])

	#Objective
	@objective(lp, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderArcSet[i]) for i in orders) + sum(c[a]*(y[a]) for a in A_hasdriver) + sum(u[a]*(w[a]) for a in A_space) )

	#Order constraints
	@constraint(lp, orderFlowBalance[i = orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))], sum(x[i,a] for a in A_minus_i[i,n]) - sum(x[i,a] for a in A_plus_i[i,n]) == 0)
	@constraint(lp, arriveDestin[i = orders], sum(sum(x[i,a] for a in A_minus_i[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in Destination[i]) == 1)
	@constraint(lp, departOrigin[i = orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), A_plus_i[i,n])) for n in Origin[i]) == 1)
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
	@constraint(lp, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in A_minus_i[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(lp, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_plus_hd[n]) == m_0[n])
	@constraint(lp, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) >= m_end[n])
	@constraint(lp, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n])- sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n],dummyarc)) for i in orders) - sum(y[a] for a in A_plus_hd[n]) == 0)

	#Linking constraints
	@constraint(lp, driverAvailability[a in A_space], sum(sum(sum(z[l,s,f] for f in fragmentscontaining[l,s,a]) for s in 1:numeffshifts) for l in 1:numlocs) == w[a] )
	for i in orders, a in orderArcSet_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in A_hasdriver_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(lp, driverStartingLocs[l in 1:numlocs, s in 1:numeffshifts], sum(sum(z[l,s,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == length(driversets[l,s]))
	@constraint(lp, driverFlowBalance[l in 1:numlocs, s in 1:numeffshifts, n in N_flow_ls[l,s]], sum(z[l,s,f] for f in F_minus_ls[l,s,n]) - sum(z[l,s,f] for f in F_plus_ls[l,s,n]) == 0)

	optimize!(lp)

	lp_obj = objective_value(lp)
	println("Fragment LP objective = ", lp_obj)
	println("Time = ", solve_time(lp))

	#bhead = getBasisHead(lp)

	return lp_obj, value.(z), value.(x), solve_time(lp)

end

#----------------------------------------------------------------------------------------#

function solvefragmentip(opt_gap, orderArcSet, orderArcSet_space, A_plus_i, A_minus_i, numeffshifts)

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60*40)
	set_optimizer_attribute(ip, "OutputFlag", 1)
	set_optimizer_attribute(ip, "MIPGap", opt_gap)

	#Variables
	@variable(ip, x[i in orders, a in orderArcSet[i]], Bin)
	@variable(ip, y[A_hasdriver] >= 0, Int)
	@variable(ip, z[l = 1:numlocs, s = 1:numeffshifts, f = 1:numfragments[l,s]] >= 0, Int)
	@variable(ip, w[a in A_space] >= 0)
	@variable(ip, ordtime[orders])

	#Objective
	@objective(ip, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderArcSet[i]) for i in orders) + sum(c[a]*(y[a] ) for a in A_hasdriver) + sum(u[a]*(w[a] ) for a in A_space) )

	#Order constraints
	@constraint(ip, orderFlowBalance[i = orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))], sum(x[i,a] for a in A_minus_i[i,n]) - sum(x[i,a] for a in A_plus_i[i,n]) == 0)
	@constraint(ip, arriveDestin[i = orders], sum(sum(x[i,a] for a in A_minus_i[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in Destination[i]) == 1)
	@constraint(ip, departOrigin[i = orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), A_plus_i[i,n])) for n in Origin[i]) == 1)
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
	@constraint(ip, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in A_minus_i[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(ip, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_plus_hd[n]) == m_0[n])
	@constraint(ip, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) >= m_end[n])
	#@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) - sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n],dummyarc)) for i in orders) - sum(y[a] for a in A_plus_hd[n]) == 0)
	@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) - sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n],dummyarc)) for i in orders) - sum(y[a] for a in A_plus_hd[n]) == 0)


	#Linking constraints
	@constraint(ip, driverAvailability[a in A_space], sum(sum(sum(z[l,s,f] for f in fragmentscontaining[l,s,a]) for s in 1:numeffshifts) for l in 1:numlocs) == w[a]  )
	for i in orders, a in orderArcSet_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in A_hasdriver_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(ip, driverStartingLocs[l in 1:numlocs, s in 1:numeffshifts], sum(sum(z[l,s,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == length(driversets[l,s]))
	@constraint(ip, driverFlowBalance[l in 1:numlocs, s in 1:numeffshifts, n in N_flow_ls[l,s]], sum(z[l,s,f] for f in F_minus_ls[l,s,n]) - sum(z[l,s,f] for f in F_plus_ls[l,s,n]) == 0)
 
	optimize!(ip)

	ip_obj = objective_value(ip)
	println("IP objective = ", ip_obj)
	println("Time = ", solve_time(ip))

	totalmiles = sum(sum(c[a]*value(x[i,a]) for a in orderArcSet[i]) for i in orders) + sum(c[a]*(value(y[a])) for a in A_hasdriver) + sum(u[a]*(value(w[a]) ) for a in A_space) 
	totaldelay = sum((value(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders)

	return ip_obj, value.(z), solve_time(ip), objective_bound(ip), value.(x), totalmiles, totaldelay
	
end

