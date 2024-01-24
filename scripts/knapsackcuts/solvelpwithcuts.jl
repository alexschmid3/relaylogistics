
function solvelpwithcuts(opt_gap, orderarcs, startcuttype)

    if startcuttype == 2
        cuttype = 1
    else
        cuttype = startcuttype
    end

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

	#Objective
	@objective(lp, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) 
		+ sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in orders) + sum(c[a]*(y[a]) for a in hasdriverarcs.A) + sum(u[a]*(w[a]) for a in primaryarcs.A_space) 
		+ lambda2 * maxhours)

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
    
    #Simple cuts
    #@constraint(lp, simplecuts[d in drivers, l = driverHomeLocs[d], s = drivershift[d], jl = tstep:tstep:24],
    # sum(z[d,f] for f in 1:numfragments[l,s] if fragworkinghours[l,s,f] == jl) <= floor(maxweeklydriverhours/jl))

    #Begin loop of cuts
    mastercuts = (vars=Dict(), rhs=Dict(), coeff=Dict())
    cutiter = 1
    originallpopt = 0
    while 1==1

        optimize!(lp)
        if cutiter == 1
            originallpopt += objective_value(lp)
        end

        println("-------------- ITER $cutiter --------------")
        println("Objective = ", objective_value(lp))

        #Knapsack cuts
        cuts = findknapsackcuts(z, cuttype)
        #cuts_nolift = findknapsackcuts(z, 4)
        #cuts_lift = findknapsackcuts(z, 5)
        
        if (cuts.numcuts == 0) & (cuttype == startcuttype) 
            println("No cuts to add!")
            println("------------------------------------")
            break
        elseif (cuts.numcuts == 0) & (cuttype == 1) & (startcuttype == 2) 
            cuttype += 1
            cuts = findknapsackcuts(z, cuttype)
            if (cuts.numcuts == 0) 
                println("No cuts to add!")
                println("------------------------------------")
                break
            end
        end

        if cuttype == 10
            varsforviz = Dict()
            for d in drivers
                varsforviz[d] = []
            end
            for i in 1:cuts.numcuts
                a,b,c = cuts.vars[i]
                push!(varsforviz[a[1]], (a[2], b[2], c[2]))
            end

            for d in drivers
                lparcs, cutarcs = [], []
                for f in 1:numfragments[driverHomeLocs[d], drivershift[d]]
                    if value(z[d,f]) > 1e-4
                        for a in 1:numarcs
                            if f in fragmentscontaining[driverHomeLocs[d],drivershift[d],a]
                                push!(lparcs, a)
                            end
                        end
                    end
                end
                for (f1,f2,f3) in varsforviz[d]
                    for f in [f1,f2,f3]
                        for a in 1:numarcs
                            if f in fragmentscontaining[driverHomeLocs[d],drivershift[d],a]
                                push!(cutarcs, a)
                            end
                        end
                    end
                end
                
                arclistlist = [driverarcs.A[d], lparcs, cutarcs]
                colorlist = [(200,200,200),(0,0,0), (255,0,0)] 
                thicknesslist = [5,11,5]
                timespacenetwork(string("outputs/viz/driver", d,"_iter", cutiter,".png"), arclistlist, colorlist, thicknesslist, 2400, 1800)
            end
        end
        
        @constraint(lp, [i in 1:cuts.numcuts], sum(cuts.coeff[i][d,j] * z[d,j] for (d,j) in cuts.vars[i]) <= cuts.rhs[i])
        cutindex = length(mastercuts.rhs)+1
        for i in 1:cuts.numcuts
            mastercuts.vars[cutindex] = cuts.vars[i]
            mastercuts.rhs[cutindex] = cuts.rhs[i]
            mastercuts.coeff[cutindex] = cuts.coeff[i]
            cutindex += 1
        end
        println("Added ", cuts.numcuts, " cuts of type $cuttype")        
        cutiter += 1
        
    end

	lp_obj = objective_value(lp)
	println("LP w/ cuts objective = ", lp_obj)
    println("Increased bound by ", round(100*(lp_obj - originallpopt)/originallpopt, digits=3), "%")
	println("Time = ", solve_time(lp))

    return lp_obj, value.(x), value.(z), solve_time(lp), objective_bound(lp), mastercuts

end