
function solvejourneymodel_paths(lprelax_flag, opt_gap, paths, delta, numeffshifts)

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60*40)
	set_optimizer_attribute(ip, "OutputFlag", 1)
	set_optimizer_attribute(ip, "MIPGap", opt_gap)

	#Variables
    if lprelax_flag == 0
        @variable(ip, x[i in orders, p in paths[i]], Bin)
        @variable(ip, y[hasdriverarcs.A] >= 0, Int)
        @variable(ip, z[d = drivers, f = 1:numfragments[driverHomeLocs[d], drivershift[d]]] >= 0, Bin)	
    elseif lprelax_flag == 1
        @variable(ip, 0 <= x[i in orders, p in paths[i]] <= 1)
	    @variable(ip, y[hasdriverarcs.A] >= 0)
	    @variable(ip, 0 <= z[d = drivers, f = 1:numfragments[driverHomeLocs[d], drivershift[d]]] <= 1)
    end
    @variable(ip, w[a in primaryarcs.A_space] >= 0)
	@variable(ip, ordtime[orders])
    @variable(ip, maxhours)

	#Objective
	@objective(ip, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) 
        + sum(sum(sum(c[a]*delta[i,a,p]*x[i,p] for a in orderarcs.A[i] ) for p in paths[i]) for i in orders) + sum(c[a]*(y[a] ) for a in hasdriverarcs.A) + sum(u[a]*(w[a] ) for a in primaryarcs.A_space) 
        + lambda2 * maxhours)

	#Order constraints
	@constraint(ip, orderpath[i in orders], sum(x[i,p] for p in paths[i]) == 1)

	#Order delivery constraints
	@constraint(ip, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in orderarcs.A_minus[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(ip, initialTrucks[n in N_0], sum(sum(sum(delta[i,a,p] * x[i,p] for a in orderarcs.A_plus[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_plus[n]) == m_0[n])
	@constraint(ip, finalTrucks[n in N_end], sum(sum(sum(delta[i,a,p] * x[i,p] for a in orderarcs.A_minus[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) >= m_end[n])
	@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(sum(delta[i,a,p] * x[i,p] for a in orderarcs.A_minus[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n],dummyarc)) for i in orders) - sum(y[a] for a in hasdriverarcs.A_plus[n]) == 0)

	#Linking constraints
	@constraint(ip, driverAvailability[a in primaryarcs.A_space], sum(sum(z[d,f] for f in fragmentscontaining[driverHomeLocs[d], drivershift[d],a]) for d in drivers) == w[a]  )
	for i in orders, p in paths[i], a in orderarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,p], -1 * delta[i,a,p])
	end
	for a in hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(ip, driverStartingLocs[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(sum(z[d,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == 1)
	@constraint(ip, driverFlowBalance[d in drivers, l = driverHomeLocs[d], s = drivershift[d], n in N_flow_ls[l,s]], sum(z[d,f] for f in F_minus_ls[l,s,n]) - sum(z[d,f] for f in F_plus_ls[l,s,n]) == 0)
	@constraint(ip, driverMaxHours[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(fragworkinghours[l,s,f] * z[d,f] for f in 1:numfragments[l,s]) <= maxhours)
	@constraint(ip, maxhours <= maxweeklydriverhours)

	#Solve
    optimize!(ip)

	#Get objective
    ip_obj = objective_value(ip)
	if lprelax_flag == 1
		println("LP objective = ", ip_obj)
	else
		println("IP objective = ", ip_obj)
	end
	println("Time = ", solve_time(ip))

	return ip_obj, value.(x), value.(z), solve_time(ip), objective_bound(ip) 
	
end


