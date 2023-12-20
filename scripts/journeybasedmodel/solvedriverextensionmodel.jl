
function solvedriverextensionmodel(lprelax_flag, opt_gap, orderArcSet, orderArcSet_space, A_plus_i, A_minus_i, numeffshifts)

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60*40)
	set_optimizer_attribute(ip, "OutputFlag", 1)
	set_optimizer_attribute(ip, "MIPGap", opt_gap)

	#Variables
    if lprelax_flag == 0
        @variable(ip, x[i in orders, a in orderArcSet[i]], Bin)
        @variable(ip, y[A_hasdriver] >= 0, Int)
        @variable(ip, z[d = drivers, f = 1:numfragments[driverHomeLocs[d], drivershift[d]]] >= 0, Int)	
    elseif lprelax_flag == 1
        @variable(ip, 0 <= x[i in orders, a in orderArcSet[i]] <= 1)
	    @variable(ip, y[A_hasdriver] >= 0)
	    @variable(ip, z[d = drivers, f = 1:numfragments[driverHomeLocs[d], drivershift[d]]] >= 0)
    end
    @variable(ip, w[a in A_space] >= 0)
	@variable(ip, ordtime[orders])
	@variable(ip, maxhours <= maxweeklydriverhours)

	#Objective
	@objective(ip, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) 
	+ sum(sum(c[a]*x[i,a] for a in orderArcSet[i]) for i in orders) + sum(c[a]*(y[a] ) for a in A_hasdriver) + sum(u[a]*w[a] for a in A_space) 
	+ lambda2 * maxhours)

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
	@constraint(ip, driverAvailability[a in A_space], sum(sum(z[d,f] for f in fragmentscontaining[driverHomeLocs[d], drivershift[d],a]) for d in drivers) == w[a]  )
	for i in orders, a in orderArcSet_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in A_hasdriver_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(ip, driverStartingLocs[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(sum(z[d,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == 1)
	@constraint(ip, driverFlowBalance[d in drivers, l = driverHomeLocs[d], s = drivershift[d], n in N_flow_ls[l,s]], sum(z[d,f] for f in F_minus_ls[l,s,n]) - sum(z[d,f] for f in F_plus_ls[l,s,n]) == 0)
	@constraint(ip, driverMaxHours[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(fragworkinghours[l,s,f] * z[d,f] for f in 1:numfragments[l,s]) <= maxhours)

	optimize!(ip)

	ip_obj = objective_value(ip)
	println("IP objective = ", ip_obj)
	println("Time = ", solve_time(ip))

	totalmiles = sum(sum(c[a]*value(x[i,a]) for a in orderArcSet[i]) for i in orders) + sum(c[a]*(value(y[a])) for a in A_hasdriver) + sum(u[a]*(value(w[a]) ) for a in A_space) 
	totaldelay = sum((value(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders)

	return ip_obj, value.(z), solve_time(ip), objective_bound(ip), value.(x), totalmiles, totaldelay
	
end

