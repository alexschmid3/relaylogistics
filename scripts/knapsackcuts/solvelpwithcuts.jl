
function solvelpwithcuts(opt_gap, orderarcs, cuttype)

    lp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(lp, "TimeLimit", 60*60*40)
	set_optimizer_attribute(lp, "OutputFlag", 0)
	set_optimizer_attribute(lp, "MIPGap", opt_gap)

	#Variables
    @variable(lp, 0 <= x[i in orders, a in orderarcs.A[i]] <= 1)
	@variable(lp, y[hasdriverarcs.A] >= 0)
	@variable(lp, 0 <= z[d = drivers, f = 1:numfragments[driverHomeLocs[d], drivershift[d]]] <= 1)
    @variable(lp, w[a in primaryarcs.A_space] >= 0)
	@variable(lp, ordtime[orders])
	@variable(lp, maxhours)
	@variable(lp, orderdelay[orders] >= 0)    #only used when deadlines turned on

	#Objective
	if deadlines_flag == 0
		@objective(lp, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in orders) + sum(c[a]*(y[a]) for a in hasdriverarcs.A) + sum(u[a]*(w[a]) for a in primaryarcs.A_space)  + lambda2 * maxhours)
    elseif deadlines_flag == 1
		@objective(lp, Min, lambda * sum(orderdelay[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in orders) + sum(c[a]*(y[a]) for a in hasdriverarcs.A) + sum(u[a]*(w[a]) for a in primaryarcs.A_space)  + lambda2 * maxhours)
		@constraint(lp, absolutedelay[i in orders], orderdelay[i] >= ordtime[i] - orderdeadline[i])
    end

	#Order constraints
	@constraint(lp, orderFlowBalance[i = orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))], sum(x[i,a] for a in orderarcs.A_minus[i,n]) - sum(x[i,a] for a in orderarcs.A_plus[i,n]) == 0)
	@constraint(lp, arriveDestin[i = orders], sum(sum(x[i,a] for a in orderarcs.A_minus[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in Destination[i]) == 1)
	@constraint(lp, departOrigin[i = orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), orderarcs.A_plus[i,n])) for n in Origin[i]) == 1)
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
	@constraint(lp, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in orderarcs.A_minus[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(lp, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_plus[n]) == m_0[n])
	@constraint(lp, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) >= m_end[n])
	@constraint(lp, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n],dummyarc)) for i in orders) - sum(y[a] for a in hasdriverarcs.A_plus[n]) == 0)

	#Linking constraints
	@constraint(lp, driverAvailability[a in primaryarcs.A_space], sum(sum(z[d,f] for f in fragmentscontaining[driverHomeLocs[d], drivershift[d],a]) for d in drivers) == w[a]  )
	for i in orders, a in orderarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(lp, driverStartingLocs[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(sum(z[d,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == 1)
	@constraint(lp, driverFlowBalance[d in drivers, l = driverHomeLocs[d], s = drivershift[d], n in N_flow_ls[l,s]], sum(z[d,f] for f in F_minus_ls[l,s,n]) - sum(z[d,f] for f in F_plus_ls[l,s,n]) == 0)
	@constraint(lp, driverMaxHours[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(fragworkinghours[l,s,f] * z[d,f] for f in 1:numfragments[l,s]) <= maxhours)
	@constraint(lp, maxhours <= maxweeklydriverhours)
    
    #Initialize cutting plane algorithm
    mastercuts = (vars=Dict(), rhs=Dict(), coeff=Dict())
    cutiter = 1
    originallpopt = 0
	starttime = time()

    #Iteratively add cuts until convergence
    while 1==1

        println("-------------- ITER $cutiter --------------")
        
        #Solve and save original objective, pre-cuts
        optimize!(lp)
        if cutiter == 1
            originallpopt += objective_value(lp)
        end
		println("Objective = ", objective_value(lp))

        #Find knapsack cuts
        cuts = findknapsackcuts(z, cuttype)
        
        #Termination condition
        if (cuts.numcuts == 0) 
            println("No cuts to add!")
            break
        end

        #Add the cuts to the LP        
        @constraint(lp, [i in 1:cuts.numcuts], sum(cuts.coeff[i][d,j] * z[d,j] for (d,j) in cuts.vars[i]) <= cuts.rhs[i])
        
        #Add the cuts to the master list
        cutindex = length(mastercuts.rhs)+1
        for i in 1:cuts.numcuts
            mastercuts.vars[cutindex] = cuts.vars[i]
            mastercuts.rhs[cutindex] = cuts.rhs[i]
            mastercuts.coeff[cutindex] = cuts.coeff[i]
            cutindex += 1
        end
        println("Added ", cuts.numcuts, " cuts of type $cuttype") 
		println("Cumulative time = ", time() - starttime, " seconds")
        
        #Iterate
        cutiter += 1
        
    end

	totalsolvetime = time() - starttime
	lp_obj = objective_value(lp)
    println("LP objective = ", originallpopt)
	println("LP w/ cuts objective = ", lp_obj)
    println("Increased bound by ", round(100*(lp_obj - originallpopt)/originallpopt, digits=3), "%")
	println("Time = ", totalsolvetime)

    #Find the LP basis arcs
	orderArcSet_basis, orderArcSet_space_basis, A_plus_i_basis, A_minus_i_basis = getnonzeroarcs(value.(x), orderarcs)
	basisarcs = (A=orderArcSet_basis, A_space=orderArcSet_space_basis, A_minus=A_minus_i_basis, A_plus=A_plus_i_basis, available=[], closelocs=[]);

    return lp_obj, value.(x), value.(z), totalsolvetime, objective_bound(lp), mastercuts, basisarcs

end