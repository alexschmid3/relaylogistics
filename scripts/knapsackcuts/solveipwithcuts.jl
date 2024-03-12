
function solveipwithcuts(opt_gap, orderarcs, numeffshifts)

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60*40)
	set_optimizer_attribute(ip, "OutputFlag", 1)
	set_optimizer_attribute(ip, "MIPGap", opt_gap)

    #Define callback function for adding knapsack cuts
    function my_callback_function(cb_data)
        
		#Get z-var values
		z_val = Dict()
        for d in drivers, f in 1:numfragments[driverHomeLocs[d], drivershift[d]]
            z_val[d,f] = callback_value(cb_data, z[d,f])
        end

		#Find cuts
        numcuts, cutvars, cutrhs, cutcoeff = findminimalcovercuts(z_val, 2)
        
		#Build constraints for IP
		for i in 1:numcuts
            con = @build_constraint( sum(cutcoeff[i][d,j] * z[d,j] for (d,j) in cutvars[i]) <= cutrhs[i] )
            MOI.submit(ip, MOI.UserCut(cb_data), con)
        end
    end

    #Set callback function
    set_attribute(ip, MOI.UserCutCallback(), my_callback_function)

	#Variables
    @variable(ip, x[i in orders, a in orderarcs.A[i]], Bin)
    @variable(ip, y[hasdriverarcs.A] >= 0, Int)
    @variable(ip, z[d = drivers, f = 1:numfragments[driverHomeLocs[d], drivershift[d]]] >= 0, Bin)
    @variable(ip, w[a in primaryarcs.A_space] >= 0)
	@variable(ip, ordtime[orders])
	@variable(ip, maxhours>=0)
	@variable(ip, orderdelay[orders] >= 0)    #only used when deadlines turned on
	
	#Objective
	if deadlines_flag == 0
		@objective(ip, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in orders) + sum(c[a]*(y[a]) for a in hasdriverarcs.A) + sum(u[a]*(w[a]) for a in primaryarcs.A_space)  + lambda2 * maxhours)
    elseif deadlines_flag == 1
		@objective(ip, Min, lambda * sum(orderdelay[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in orders) + sum(c[a]*(y[a]) for a in hasdriverarcs.A) + sum(u[a]*(w[a]) for a in primaryarcs.A_space)  + lambda2 * maxhours)
		@constraint(ip, absolutedelay[i in orders], orderdelay[i] >= ordtime[i] - orderdeadline[i])
    end

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

	#Linking constraints
	@constraint(ip, driverAvailability[a in primaryarcs.A_space], sum(sum(z[d,f] for f in fragmentscontaining[driverHomeLocs[d], drivershift[d],a]) for d in drivers) == w[a]  )
	for i in orders, a in orderarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(ip, driverStartingLocs[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(sum(z[d,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == 1)
	@constraint(ip, driverFlowBalance[d in drivers, l = driverHomeLocs[d], s = drivershift[d], n in N_flow_ls[l,s]], sum(z[d,f] for f in F_minus_ls[l,s,n]) - sum(z[d,f] for f in F_plus_ls[l,s,n]) == 0)
	@constraint(ip, driverMaxHours[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(fragworkinghours[l,s,f] * z[d,f] for f in 1:numfragments[l,s]) <= maxhours)
	@constraint(ip, maxhours <= maxweeklydriverhours)

	#Symmetry breaking constraints
	if symmetrybreaking_flag == 1
		for l in 1:numlocs, s in 1:numeffshifts
			for j in 1:length(driversets[l,s])-1
				d1, d2 = driversets[l,s][j], driversets[l,s][j+1]
				@constraint(ip, sum(z[d1,f] for f in setdiff(1:numfragments[l,s], workingfragments[l,s])) <= sum(z[d2,f] for f in setdiff(1:numfragments[l,s], workingfragments[l,s])) )
			end
		end
	end

	optimize!(ip)

	ip_obj = objective_value(ip)
	println("IP objective = ", ip_obj)
	println("Time = ", solve_time(ip))
    
	return ip_obj, value.(x), value.(z), solve_time(ip), objective_bound(ip)
	
end

