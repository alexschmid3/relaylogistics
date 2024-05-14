
function solvejourneymodel(lprelax_flag, opt_gap, arcspassed, currentdatetime)
			
	totaldelta = Dates.value(Dates.Hour(currentdatetime - weekstart))

	if arcspassed == -1
		journeysfor = Dict()
		for (i1,i2,i3,i4) in currfragments.driversets
			journeysfor[i1,i2,i3,i4] = [j for j in 1:currfragments.numfragments[i1,i2,i3,i4]]
		end
	else
		journeysfor = Dict()
		for (i1,i2,i3,i4) in currfragments.driversets
			journeysfor[i1,i2,i3,i4] = arcspassed[i1,i2,i3,i4]
		end
	end

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60)
	set_optimizer_attribute(ip, "OutputFlag", 0)
	set_optimizer_attribute(ip, "MIPGap", opt_gap)

	A_space_all = primaryarcs.A_space
	for i in currstate.orders
		goodones = [a for a in currarcs.orderarcs.A_space[i] if nodesLookup[arcLookup[a][1]][2] < horizon]
		A_space_all = union(A_space_all, goodones)
	end

	#Variables
	if lprelax_flag == 0
		@variable(ip, x[i in currstate.orders, a in currarcs.orderarcs.A[i]], Bin)
		@variable(ip, y[currarcs.hasdriverarcs.A] >= 0, Int)
		@variable(ip, z[(i1,i2,i3,i4) = currfragments.driversets, f = journeysfor[i1,i2,i3,i4]] >= 0, Int)	
	elseif lprelax_flag == 1
		@variable(ip, 0 <= x[i in currstate.orders, a in currarcs.orderarcs.A[i]] <= 1)
		@variable(ip, y[currarcs.hasdriverarcs.A] >= 0)
		@variable(ip, z[(i1,i2,i3,i4) = currfragments.driversets, f = journeysfor[i1,i2,i3,i4]] >= 0)
	end
	@variable(ip, w[a in A_space_all] >= 0)
	@variable(ip, ordtime[currstate.orders])

	#Objective
	@objective(ip, Min, lambda * sum(ordtime[i] for i in currstate.orders) #sum((ordtime[i] - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i] for i in currstate.orders) 
	+ sum(sum(c[a]*x[i,a] for a in currarcs.orderarcs.A[i]) for i in currstate.orders) + sum(c[a]*(y[a]) for a in currarcs.hasdriverarcs.A) + sum(u[a]*(w[a] ) for a in primaryarcs.A_space) )
	#@objective(ip, Min, lambda * sum(ordtime[i] for i in [2]) #sum((ordtime[i] - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i] for i in currstate.orders) 
	#+ sum(sum(c[a]*x[i,a] for a in currarcs.orderarcs.A[i]) for i in [2]) ) #+ sum(c[a]*(y[a]) for a in currarcs.hasdriverarcs.A) + sum(u[a]*(w[a] ) for a in primaryarcs.A_space) )

	#Order constraints
	@constraint(ip, orderFlowBalance[i = currstate.orders, n in setdiff([n2 for n2 in 1:numnodes], union(currstate.Origin[i], currstate.Destination[i]))], sum(x[i,a] for a in currarcs.orderarcs.A_minus[i,n]) - sum(x[i,a] for a in currarcs.orderarcs.A_plus[i,n]) == 0)
	@constraint(ip, arriveDestin[i = currstate.orders], sum(sum(x[i,a] for a in currarcs.orderarcs.A_minus[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in currstate.Destination[i]) == 1)
	@constraint(ip, departOrigin[i = currstate.orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), currarcs.orderarcs.A_plus[i,n])) for n in currstate.Origin[i]) == 1)
	for i in setdiff(currstate.orders, currstate.ordersinprogress)
		extendedorderarc = extendedarcs[last(currstate.Origin[i]), last(currstate.Destination[i])]
		if extendedorderarc in currarcs.orderarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end

	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(currstate.orders, currstate.ordersinprogress), n in currstate.Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
		if a in currarcs.orderarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(ip, deliveryTime[i in currstate.orders], ordtime[i] - sum(sum((totaldelta + arcfinishtime[a]) * x[i,a] for a in currarcs.orderarcs.A_minus[i,n]) for n in currstate.Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(ip, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(currarcs.orderarcs.A_plus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == currstate.m_0[n])
	@constraint(ip, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(currarcs.orderarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) >= currstate.m_end[n])
	#@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(currarcs.orderarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(currarcs.orderarcs.A_plus[i,n],dummyarc)) for i in currstate.orders) - sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == 0)
	@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(currarcs.orderarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(currarcs.orderarcs.A_plus[i,n],dummyarc)) for i in currstate.orders) - sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == 0)

	#Trucks in transit --> hold the intransit trucks at their current loc until they "become available"
	for item in currstate.trucksintransit
		l, availabletime, numtrucksheld = item
		for t in 0:tstep:availabletime - tstep
			@constraint(ip, y[arcs[nodes[l,t], nodes[l,t+tstep]]] >= min(currstate.m_0[l],numtrucksheld))
		end
	end

	#Linking constraints
	@constraint(ip, driverAvailability[a in A_space_all], sum(sum(z[(i1,i2,i3,i4),f] for f in intersect(currfragments.fragmentscontaining[i1,i2,i3,i4,a], journeysfor[i1,i2,i3,i4]); init=0) for (i1,i2,i3,i4) in currfragments.driversets; init=0) == w[a]  )
	for i in currstate.orders, a in [a for a in currarcs.orderarcs.A_space[i] if nodesLookup[arcLookup[a][1]][2] < horizon]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in currarcs.hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(ip, driverStartingLocs[(i1,i2,i3,i4) in currfragments.driversets], sum(sum(z[(i1,i2,i3,i4),f] for f in intersect(currfragments.F_plus_g[i1,i2,i3,i4,n], journeysfor[i1,i2,i3,i4])) for n in [i3]) == length(currfragments.driversingroup[i1,i2,i3,i4]))
	@constraint(ip, driverFlowBalance[(i1,i2,i3,i4) in currfragments.driversets, n in currfragments.N_flow_g[i1,i2,i3,i4]], sum(z[(i1,i2,i3,i4),f] for f in intersect(currfragments.F_minus_g[i1,i2,i3,i4,n], journeysfor[i1,i2,i3,i4])) - sum(z[(i1,i2,i3,i4),f] for f in intersect(currfragments.F_plus_g[i1,i2,i3,i4,n], journeysfor[i1,i2,i3,i4])) == 0)

	#--------------------------------------------#

	optimize!(ip)

	#--------------------------------------------#

	ip_obj = objective_value(ip)
	if lprelax_flag == 1
		println("LP objective = ", ip_obj)
	else
		println("IP objective = ", ip_obj)
	end
	println("Time = ", solve_time(ip))

	totalmiles = sum(sum(c[a]*value(x[i,a]) for a in currarcs.orderarcs.A[i]) for i in currstate.orders) + sum(c[a]*(value(y[a])) for a in currarcs.hasdriverarcs.A) + sum(u[a]*(value(w[a]) ) for a in primaryarcs.A_space) 
	totaldelay = sum((value(ordtime[i]) - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i] for i in currstate.orders)

	#--------------------------------------------#

	useddriverjourneys = Dict()
	for (i1,i2,i3,i4) in currfragments.driversets
		useddriverjourneys[i1,i2,i3,i4] = []
	end
	for (i1,i2,i3,i4) in currfragments.driversets, j in journeysfor[i1,i2,i3,i4]
		if value(z[(i1,i2,i3,i4),j]) > 1e-4
			push!(useddriverjourneys[i1,i2,i3,i4], j)
		end
	end

	#--------------------------------------------#

	return ip_obj, value.(x), value.(z), value.(w), value.(y), solve_time(ip), objective_bound(ip), useddriverjourneys #, totalmiles, totaldelay
	
end

#=
myarcs, usedarcs = [], []
for j in 1:currfragments.numfragments[14, 2, 23, -24]
	myarcs = union(myarcs, currfragments.fragmentarcs[14, 2, 23, -24, j])
end
for j in 1:currfragments.numfragments[14, 2, 23, -24]
	if value(z[(14, 2, 23, -24), j]) > 1e-4
		usedarcs = union(usedarcs, currfragments.fragmentarcs[14, 2, 23, -24, j])
	end
end
timespacenetwork("outputs/viz/aaa_all.png", [myarcs, usedarcs], [(150,150,150),(0,0,0)], [3,6,6], ["solid","solid","solid"], [0,0,0], 2400, 1800)

=#
