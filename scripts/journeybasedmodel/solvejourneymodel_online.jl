
function solvejourneymodel(lprelax_flag, opt_gap, currstate, currarcs, currfragments, orderarcs)

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60*40)
	set_optimizer_attribute(ip, "OutputFlag", 1)
	set_optimizer_attribute(ip, "MIPGap", opt_gap)

	#Variables
    if lprelax_flag == 0
        @variable(ip, x[i in currstate.orders, a in orderarcs.A[i]], Bin)
        @variable(ip, y[currarcs.hasdriverarcs.A] >= 0, Int)
        @variable(ip, z[l = 1:numlocs, s = 1:numeffshifts, f = 1:currfragments.numfragments[l,s]] >= 0, Int)	
    elseif lprelax_flag == 1
        @variable(ip, 0 <= x[i in currstate.orders, a in orderarcs.A[i]] <= 1)
	    @variable(ip, y[currarcs.hasdriverarcs.A] >= 0)
	    @variable(ip, z[l = 1:numlocs, s = 1:numeffshifts, f = 1:currfragments.numfragments[l,s]] >= 0)
    end
    @variable(ip, w[a in primaryarcs.A_space] >= 0)
	@variable(ip, ordtime[currstate.orders])

	#Objective
	@objective(ip, Min, lambda * sum((ordtime[i] - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i] for i in currstate.orders) 
	+ sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in currstate.orders) + sum(c[a]*(y[a] ) for a in currarcs.hasdriverarcs.A) + sum(u[a]*(w[a] ) for a in primaryarcs.A_space) )

	#Order constraints
	@constraint(ip, orderFlowBalance[i = currstate.orders, n in setdiff([n2 for n2 in 1:numnodes], union(currstate.Origin[i], currstate.Destination[i]))], sum(x[i,a] for a in orderarcs.A_minus[i,n]) - sum(x[i,a] for a in orderarcs.A_plus[i,n]) == 0)
	@constraint(ip, arriveDestin[i = currstate.orders], sum(sum(x[i,a] for a in orderarcs.A_minus[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in currstate.Destination[i]) == 1)
	@constraint(ip, departOrigin[i = currstate.orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), orderarcs.A_plus[i,n])) for n in currstate.Origin[i]) == 1)
	for i in setdiff(currstate.orders, currstate.ordersinprogress)
		extendedorderarc = extendedarcs[last(currstate.Origin[i]), last(currstate.Destination[i])]
		if extendedorderarc in orderarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end
	
	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(currstate.orders, currstate.ordersinprogress), n in currstate.Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
		if a in orderarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(ip, deliveryTime[i in currstate.orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in orderarcs.A_minus[i,n]) for n in currstate.Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(ip, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == currstate.m_0[n])
	@constraint(ip, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) >= currstate.m_end[n])
	#@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n],dummyarc)) for i in currstate.orders) - sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == 0)
	@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n],dummyarc)) for i in currstate.orders) - sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == 0)

	#Linking constraints
	@constraint(ip, driverAvailability[a in primaryarcs.A_space], sum(sum(sum(z[l,s,f] for f in currfragments.fragmentscontaining[l,s,a]) for s in 1:numshifts) for l in 1:numlocs) == w[a]  )
	for i in currstate.orders, a in orderarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in currarcs.hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(ip, driverStartingLocs[l in 1:numlocs, s in 1:numeffshifts], sum(sum(z[l,s,f] for f in currfragments.F_plus_ls[l,s,n]) for n in currfragments.driverSetStartNodes[l,s]) == length(currfragments.driversets[l,s]))
	@constraint(ip, driverFlowBalance[l in 1:numlocs, s in 1:numeffshifts, n in currfragments.N_flow_ls[l,s]], sum(z[l,s,f] for f in currfragments.F_minus_ls[l,s,n]) - sum(z[l,s,f] for f in currfragments.F_plus_ls[l,s,n]) == 0)

	optimize!(ip)

	ip_obj = objective_value(ip)
	if lprelax_flag == 1
		println("LP objective = ", ip_obj)
	else
		println("IP objective = ", ip_obj)
	end
	println("Time = ", solve_time(ip))

	totalmiles = sum(sum(c[a]*value(x[i,a]) for a in orderarcs.A[i]) for i in currstate.orders) + sum(c[a]*(value(y[a])) for a in currarcs.hasdriverarcs.A) + sum(u[a]*(value(w[a]) ) for a in primaryarcs.A_space) 
	totaldelay = sum((value(ordtime[i]) - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i] for i in currstate.orders)

	return ip_obj, value.(x), value.(z), solve_time(ip), objective_bound(ip) #, totalmiles, totaldelay
	
end

