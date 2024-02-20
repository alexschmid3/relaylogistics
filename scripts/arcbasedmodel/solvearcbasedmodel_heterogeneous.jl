
function solvearcbasedmodel(orderarcs, lprelaxflag)

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60*47)
	set_optimizer_attribute(ip, "OutputFlag", 1)

	#Variables
	if lprelaxflag == 0
		@variable(ip, x[i in orders, a in orderarcs.A[i]], Bin)
		@variable(ip, y[hasdriverarcs.A] >= 0, Int)
		@variable(ip, z[d = drivers, driverarcs.A[d]], Bin)
		@variable(ip, w[a in primaryarcs.A_space] >= 0)
		@variable(ip, ordtime[orders])
		@variable(ip, maxhours>=0)
	elseif lprelaxflag == 1
		@variable(ip, 0 <= x[i in orders, a in orderarcs.A[i]] <= 1)
		@variable(ip, y[hasdriverarcs.A] >= 0)
		@variable(ip, 0 <= z[d = drivers, driverarcs.A[d]] <= 1)
		@variable(ip, w[a in primaryarcs.A_space] >= 0)
		@variable(ip, ordtime[orders])
		@variable(ip, maxhours>=0)
	end

	#Objective
	@objective(ip, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) 
		+ sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in orders) + sum(c[a]*y[a] for a in hasdriverarcs.A) + sum(u[a]*w[a] for a in primaryarcs.A_space) 
		+ lambda2 * maxhours)
	
	#Order constraints
	@constraint(ip, orderFlowBalance[i = orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))], sum(x[i,a] for a in orderarcs.A_minus[i,n]) - sum(x[i,a] for a in orderarcs.A_plus[i,n]) == 0)
	@constraint(ip, arriveDestin[i = orders], sum(sum(x[i,a] for a in orderarcs.A_minus[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in Destination[i]) == 1)
	@constraint(ip, departOrigin[i = orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), orderarcs.A_plus[i,n])) for n in Origin[i]) == 1)
	for i in setdiff(orders, ordersinprogress)
		extendedorderarc = extendedarcs[last(Origin[i]), last(Destination[i])]
		if extendedorderarc in orderarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end
	
	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(orders, ordersinprogress), n in Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
		if a in orderarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(ip, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in orderarcs.A_minus[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(ip, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_plus[n]) == m_0[n])
	@constraint(ip, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) >= m_end[n])
	@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n],dummyarc)) for i in orders) - sum(y[a] for a in hasdriverarcs.A_plus[n]) == 0)

	#Driver constraints
	@constraint(ip, driverAvailability[a in primaryarcs.A_space], sum(z[d,a] for d in intersect(drivers, driverarcs.available[a])) == w[a])
	for i in orders, a in orderarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end
	@constraint(ip, driverFlowBalance[d in drivers, n in N_flow_d[d]], sum(z[d,a] for a in driverarcs.A_minus[d,n]) - sum(z[d,a] for a in driverarcs.A_plus[d,n]) == 0)
	@constraint(ip, driverStartingLocs[d in drivers], sum(sum(z[d,a] for a in driverarcs.A_plus[d,n]) for n in driverStartNodes[d]) == 1)
	@constraint(ip, returnHome[d in drivers, arcset in R_off[d]], sum(z[d,a] for a in arcset) >= 1)  

	#Driver hours
	@constraint(ip, driverMaxHours[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum((nodesLookup[arcLookup[a][2]][2]-nodesLookup[arcLookup[a][1]][2]) * z[d,a] for a in driverarcs.A_space[d]) <= maxhours)
	@constraint(ip, maxhours <= maxweeklydriverhours)

	optimize!(ip)

	ip_obj = objective_value(ip)
	if lprelaxflag == 1
		println("LP objective = ", ip_obj)
	else
		println("IP objective = ", ip_obj)
	end
	println("Time = ", solve_time(ip))

	return ip_obj, value.(x), value.(z), solve_time(ip), objective_bound(ip)

end