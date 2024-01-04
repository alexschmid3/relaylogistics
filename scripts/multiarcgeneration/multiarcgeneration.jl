
include("preprocessmagsets.jl")

function sparsemasterproblem(magarcs, hasdriverarcs, timelimit)

	smp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(smp, "TimeLimit", timelimit)
	set_optimizer_attribute(smp, "OutputFlag", 0)
    set_optimizer_attribute(smp, "MIPGap", 0.001)

	#Variables
	x = Dict()
	for i in orders, a in magarcs.A[i]
	    global x[i,a] = @variable(smp, lower_bound = 0) #, upper_bound = 1) <-- UB implied by order origin/dest con & messes up reduced cost sign
	    set_name(x[i,a], string("x[",i,",",a,"]")) 
	end
	@variable(smp, y[hasdriverarcs.A] >= 0)
    @variable(smp, 0 <= z[d = drivers, f = 1:numfragments[driverHomeLocs[d], drivershift[d]]] <= 1)
	@variable(smp, w[a in primaryarcs.A_space] >= 0)
	@variable(smp, ordtime[orders])
    @variable(smp, maxhours)

	#Objective
	@objective(smp, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in magarcs.A[i]) for i in orders) + sum(c[a]*y[a] for a in hasdriverarcs.A) + sum(u[a]*w[a] for a in primaryarcs.A_space) + lambda2 * maxhours)
		
	#Order constraints
	@constraint(smp, orderFlowBalance[i = orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))], sum(x[i,a] for a in magarcs.A_minus[i,n]) - sum(x[i,a] for a in magarcs.A_plus[i,n]) == 0)
	@constraint(smp, arriveDestin[i = orders], sum(sum(x[i,a] for a in magarcs.A_minus[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in Destination[i]) == 1)
	@constraint(smp, departOrigin[i = orders], sum(sum(x[i,a] for a in intersect(union(primaryarcs.A_space, dummyarc), magarcs.A_plus[i,n])) for n in Origin[i]) == 1)
	for i in setdiff(orders, ordersinprogress)
		extendedorderarc = extendedarcs[last(Origin[i]), last(Destination[i])]
		if extendedorderarc in magarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end
	
	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(orders, ordersinprogress), n in Origin[i], a in setdiff(primaryarcs.A_plus[n], union(primaryarcs.A_space, dummyarc))
		if a in magarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(smp, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in magarcs.A_minus[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(smp, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(magarcs.A_plus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_plus[n]) == m_0[n])
	@constraint(smp, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(magarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) >= m_end[n])
	@constraint(smp, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(magarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n])- sum(sum(x[i,a] for a in setdiff(magarcs.A_plus[i,n],dummyarc)) for i in orders) - sum(y[a] for a in hasdriverarcs.A_plus[n]) == 0)

	#Linking constraints
	@constraint(smp, driverAvailability[a in primaryarcs.A_space], sum(sum(z[d,f] for f in fragmentscontaining[driverHomeLocs[d], drivershift[d],a]) for d in drivers) == w[a]  )
	for i in orders, a in magarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
    @constraint(smp, driverStartingLocs[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(sum(z[d,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == 1)
	@constraint(smp, driverFlowBalance[d in drivers, l = driverHomeLocs[d], s = drivershift[d], n in N_flow_ls[l,s]], sum(z[d,f] for f in F_minus_ls[l,s,n]) - sum(z[d,f] for f in F_plus_ls[l,s,n]) == 0)
	@constraint(smp, driverMaxHours[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(fragworkinghours[l,s,f] * z[d,f] for f in 1:numfragments[l,s]) <= maxhours)
	@constraint(smp, maxhours <= maxweeklydriverhours)

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

function converttosparsearray(jumpcontainer::JuMP.Containers.SparseAxisArray, keysize1::Int, keysize2::Int)

    if keysize2 == 1
        newarray = spzeros(keysize1)
        for indexkey in eachindex(jumpcontainer)
            newarray[indexkey] = jumpcontainer[indexkey]
        end
    else
        newarray = spzeros(keysize1, keysize2)
        for indexkey in eachindex(jumpcontainer)
            newarray[indexkey[1],indexkey[2]] = jumpcontainer[indexkey]
        end
    end

    return newarray

end

#----------------------------------------------------------------------------------------#

function getdualvalues(smpconstraints, setvariables)

    alpha = converttosparsearray(dual.(smpconstraints.con_orderFlowBalance), numorders, extendednumnodes)
    beta = Array(dual.(smpconstraints.con_departOrigin))
    gamma = Array(dual.(smpconstraints.con_arriveDestin))
    theta = Array(dual.(smpconstraints.con_initialTrucks))
    nu = Array(dual.(smpconstraints.con_finalTrucks))
    mu = Array(dual.(smpconstraints.con_truckFlowBalance))
    xi = Array(dual.(smpconstraints.con_driverAvailability))
    psi = Array(dual.(smpconstraints.con_deliveryTime))

    #if length(setvariables) > 0
    #    zeta = converttosparsearray(dual.(setvars_con2), numorders, extendednumarcs)
    #    return alpha, beta, gamma, theta, nu, mu, xi, psi, zeta
    #else 
    return alpha, beta, gamma, theta, nu, mu, xi, psi
    #end

end

#----------------------------------------------------------------------------------------#

function findarcvariablereducedcosts_original(M, smpconstraints, setvariables)

    alpha = dual.(smpconstraints.con_orderFlowBalance)
    beta = dual.(smpconstraints.con_departOrigin)
    gamma = dual.(smpconstraints.con_arriveDestin)
    theta = dual.(smpconstraints.con_initialTrucks)
    nu = dual.(smpconstraints.con_finalTrucks)
    mu = dual.(smpconstraints.con_truckFlowBalance)
    xi = dual.(smpconstraints.con_driverAvailability)
    psi = dual.(smpconstraints.con_deliveryTime)

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

#=
include("scripts/multiarcgeneration/multiarcgeneration.jl")
include("scripts/multiarcgeneration/preprocessmagsets.jl")
arcredcosts = findarcvariablereducedcosts(M, smpconstraints, setvariables)
arcredcosts_corr = findarcvariablereducedcosts_original(M, smpconstraints, setvariables)

wrongcounter = 0
for i in orders, a in setdiff(orderarcs.A[i], 3392)
    if abs(arcredcosts[i,a] - arcredcosts_corr[i,a]) > 1e-4
        println("$i, $a --> ", arcredcosts_corr[i,a], " vs. ", arcredcosts[i,a])

        wrongcounter += 1
        if wrongcounter > 100
            break
        end
    end
end

wrongcounter = 0
for i in orders, a in setdiff(magarcs.A[i], 3392)
    if abs(arcredcosts[i,a] - reduced_cost(x[i,a])) > 1e-4
        println("$i, $a --> ", reduced_cost(x[i,a]), " vs. ", arcredcosts[i,a])

        wrongcounter += 1
        if wrongcounter > 100
            break
        end
    end
end

=#

#----------------------------------------------------------------------------------------#

### TEST FOR SPEED (MAYBE SCALE UP INSTANCE?)

function findarcvariablereducedcosts(M, smpconstraints, setvariables)

    alpha, beta, gamma, theta, nu, mu, xi, psi = getdualvalues(smpconstraints, setvariables, 0)
    
    arcredcosts = zeros(numorders, extendednumarcs)
    for i in orders
        arcredcosts[i,:] += c[1:extendednumarcs] + M.alpha[i] * alpha[i,:] + M.beta[i] * beta[i] + M.gamma[i] * gamma[i] + M.psi[i] * psi[i]
        arcredcosts[i,:] += M.theta * theta + M.nu * nu + M.mu * mu + M.xi * xi
    end

	return arcredcosts

end

#----------------------------------------------------------------------------------------#

### TEST FOR SPEED (MAYBE SCALE UP INSTANCE?)

function findarcvariablereducedcosts_varsetting(M, smpconstraints, setvariables)

    alpha, beta, gamma, theta, nu, mu, xi, psi = getdualvalues(smpconstraints, setvariables)
    
    arcredcosts = zeros(numorders, extendednumarcs)
    for i in orders
        arcredcosts[i,:] += c[1:extendednumarcs] + M.alpha[i] * alpha[i,:] + M.beta[i] * beta[i] + M.gamma[i] * gamma[i] + M.psi[i] * psi[i]
        arcredcosts[i,:] += M.theta * theta + M.nu * nu + M.mu * mu + M.xi * xi
    end

    for (i,a) in setvariables
        arcredcosts[i,a] = 0.0 
    end

	return arcredcosts

end

#----------------------------------------------------------------------------------------#

function sagsubproblem(i, arcredcosts)

    spstarttime = time()

    mostnegativearc = setdiff(orderarcs.A[i], dummyarc)[argmin([arcredcosts[i,a] for a in setdiff(orderarcs.A[i], dummyarc)])]
    minreducedcost = arcredcosts[i, mostnegativearc]
    shortestpatharcs = [mostnegativearc]
    shortestpathnodes = []

    return minreducedcost, shortestpathnodes, shortestpatharcs, time() - spstarttime

end

#----------------------------------------------------------------------------------------#

function magsubproblem(i, arcredcosts, subproblemsets)

    spstarttime = time()

	cSP = spzeros(maximum(subproblemsets.arclist[i]))
	for a in subproblemsets.arclist[i]
		if a <= extendednumarcs
			cSP[a] = arcredcosts[i,a]
		end	
	end

	#Initialize shortest path algorithm 
	currdistance = repeat([999999999.0],outer=[subproblemsets.numnodes])
	currdistance[subproblemsets.dummyorig] = 0
	prevnode, prevarc = zeros(subproblemsets.numnodes), zeros(subproblemsets.numnodes)
	
	#Loop over time-space arcs in order of start time
	for a in subproblemsets.arclist[i]
		n_end, n_start = subproblemsets.arclookup[i][a][2], subproblemsets.arclookup[i][a][1]
		if currdistance[n_end] > currdistance[n_start] + cSP[a] + 0.000001
			currdistance[n_end] = currdistance[n_start] + cSP[a]
			prevnode[n_end] = n_start
			prevarc[n_end] = a
		end
	end

	#Format the shortest path output
	shortestpathnodes_rev = [subproblemsets.dummydest]
	shortestpatharcs_rev = []
	node = subproblemsets.dummydest
	while node != subproblemsets.dummyorig
		push!(shortestpatharcs_rev, Int(prevarc[node]))
		node = Int(prevnode[node])
		push!(shortestpathnodes_rev, node)
	end
	shortestpathnodes = reverse(shortestpathnodes_rev[2:length(shortestpathnodes_rev)-1]) 
	shortestpatharcs = reverse(shortestpatharcs_rev[2:length(shortestpatharcs_rev)-1]) 

	return currdistance[subproblemsets.dummydest], shortestpathnodes, shortestpatharcs, time() - spstarttime

end

#----------------------------------------------------------------------------------------#

function updatearcsets(magarcs, addarcs)

    newarcs = []

    for (i,a) in addarcs
        if !(a in magarcs.A[i])
            push!(magarcs.A[i], a)
            if a in primaryarcs.A_space
                push!(magarcs.A_space[i], a)
            end
            n_plus = arcLookup[a][1]
            n_minus = arcLookup[a][2]
            push!(magarcs.A_plus[i,n_plus], a)
            push!(magarcs.A_minus[i,n_minus], a)
            push!(newarcs, (i,a))
        end
    end

    return magarcs, newarcs

end

#----------------------------------------------------------------------------------------#

function multiarcgeneration!(magarcs, variablefixingthreshold, hasdriverarcs)

	#Build sparse master problem
	#smp, x, y, z, w, smpconstraints = sparsemasterproblem(magarcs, hasdriverarcs, 60*60)

    smp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(smp, "TimeLimit", 60*60)
	set_optimizer_attribute(smp, "OutputFlag", 0)
    set_optimizer_attribute(smp, "MIPGap", 0.001)

	#Variables
	x = Dict()
	for i in orders, a in magarcs.A[i]
	    global x[i,a] = @variable(smp, lower_bound = 0) #, upper_bound = 1) <-- UB implied by order origin/dest con & messes up reduced cost sign
	    set_name(x[i,a], string("x[",i,",",a,"]")) 
	end
	@variable(smp, y[hasdriverarcs.A] >= 0)
    @variable(smp, 0 <= z[d = drivers, f = 1:numfragments[driverHomeLocs[d], drivershift[d]]] <= 1)
	@variable(smp, w[a in primaryarcs.A_space] >= 0)
	@variable(smp, ordtime[orders])
    @variable(smp, maxhours)

	#Objective
	@objective(smp, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in magarcs.A[i]) for i in orders) + sum(c[a]*y[a] for a in hasdriverarcs.A) + sum(u[a]*w[a] for a in primaryarcs.A_space) + lambda2 * maxhours)
		
	#Order constraints
	@constraint(smp, orderFlowBalance[i = orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))], sum(x[i,a] for a in magarcs.A_minus[i,n]) - sum(x[i,a] for a in magarcs.A_plus[i,n]) == 0)
	@constraint(smp, arriveDestin[i = orders], sum(sum(x[i,a] for a in magarcs.A_minus[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in Destination[i]) == 1)
	@constraint(smp, departOrigin[i = orders], sum(sum(x[i,a] for a in intersect(union(primaryarcs.A_space, dummyarc), magarcs.A_plus[i,n])) for n in Origin[i]) == 1)
	for i in setdiff(orders, ordersinprogress)
		extendedorderarc = extendedarcs[last(Origin[i]), last(Destination[i])]
		if extendedorderarc in magarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end
	
	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(orders, ordersinprogress), n in Origin[i], a in setdiff(primaryarcs.A_plus[n], union(primaryarcs.A_space, dummyarc))
		if a in magarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(smp, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in magarcs.A_minus[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(smp, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(magarcs.A_plus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_plus[n]) == m_0[n])
	@constraint(smp, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(magarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) >= m_end[n])
	@constraint(smp, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(magarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n])- sum(sum(x[i,a] for a in setdiff(magarcs.A_plus[i,n],dummyarc)) for i in orders) - sum(y[a] for a in hasdriverarcs.A_plus[n]) == 0)

	#Linking constraints
	@constraint(smp, driverAvailability[a in primaryarcs.A_space], sum(sum(z[d,f] for f in fragmentscontaining[driverHomeLocs[d], drivershift[d],a]) for d in drivers) == w[a]  )
	for i in orders, a in magarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
    @constraint(smp, driverStartingLocs[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(sum(z[d,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == 1)
	@constraint(smp, driverFlowBalance[d in drivers, l = driverHomeLocs[d], s = drivershift[d], n in N_flow_ls[l,s]], sum(z[d,f] for f in F_minus_ls[l,s,n]) - sum(z[d,f] for f in F_plus_ls[l,s,n]) == 0)
	@constraint(smp, driverMaxHours[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(fragworkinghours[l,s,f] * z[d,f] for f in 1:numfragments[l,s]) <= maxhours)
	@constraint(smp, maxhours <= maxweeklydriverhours)

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

	#------------------------------------------------------#

	#Initialize column generation 
	cg_iter = 1
	smpobjectives, smptimes, pptimes, pptimes_par, lowerbounds = [], [], [], [], []
	listlength = convert(Int64, ceil(length(orders)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)

	#------------------------------------------------------#

	#Pre-processing
    M, subproblemsets = preprocessmagsets(orderarcs.A);

	#------------------------------------------------------#

    setvariables = []

    #------------------------------------------------------#

	while cg_iter <= 100000

		#-----------SOLVE SMP-----------#

		println("-------- ITERATION $cg_iter --------")

		status = optimize!(smp)
		if termination_status(smp) != MOI.OPTIMAL
			println(termination_status(smp))
			return 100000000, smp, x, y, z, w, magarcs.A
		end
		smpobj, smptime = objective_value(smp), solve_time(smp)
		
		push!(smpobjectives, copy(smpobj))
		push!(smptimes, copy(smptime))	
		println("Solved SMP with objective = ", smpobj, " in iteration $cg_iter (", sum(length(magarcs.A[i]) for i in orders), " arcs)")

        #------------SUBPROBLEMS------------#

		#Calculate reduced costs
        if setvariables == []
            arcredcosts = findarcvariablereducedcosts(M, smpconstraints, setvariables)
        else
            arcredcosts = findarcvariablereducedcosts_varsetting(M, smpconstraints, setvariables)
        end

		#Run shortest path for each order to find new arcs
		dptimelist = []
		addarcs, minreducedcosts = [], []
		for i in orders
			if onearcatatime_flag == 0	
				minreducedcost, shortestpathnodes, shortestpatharcs, sptime = magsubproblem(i, arcredcosts, subproblemsets)
			elseif onearcatatime_flag == 1
				minreducedcost, shortestpathnodes, shortestpatharcs, sptime = sagsubproblem(i, arcredcosts)
			end

			push!(dptimelist, sptime)
			push!(minreducedcosts, minreducedcost)

			if minreducedcost < -0.0001
				for a in shortestpatharcs
					push!(addarcs, (i,a))
				end				
			end
		end
		
		#"Parallelize" subproblem times
		shuffleddptimes = shuffle_partition(length(orders))
		dptimelistsums = [sum(dptimelist[j_inner] for j_inner in shuffleddptimes[j_outer]) for j_outer in 1:length(shuffleddptimes)]
		push!(pptimes_par, maximum(dptimelistsums))
        push!(pptimes, sum(dptimelist))
		
        #Update the lowerbound 
        try
			push!(lowerbounds, smpobj + sum(minreducedcosts[k] for k in 1:length(minreducedcosts) if minreducedcosts[k] <= -0.000001))
		catch
			push!(lowerbounds, smpobj)
		end

		#------COUNT ARCS AND PATHS-----#

		if saveconvergencedata_flag >= 0
			totalorderarcs = sum(length(magarcs.A[i]) for i in orders)
			totalorderpaths = sum([findallpaths(magarcs.A_plus, i) for i in orders])
			maximprove = minimum(minreducedcosts) * sum(sum(value(x[i,a]) for a in magarcs.A[i]) for i in orders)
			write_cg_conv(convergencedatafilename, cg_iter, maximprove, totalorderarcs, totalorderpaths, smpobj)
		end
	
		#-------ADD NEW VARIABLES-------#

		#Add new arcs to order arc sets
		magarcs, newarcs = updatearcsets(magarcs, addarcs)

		#Add new arcs to model as x-variables
		for (i,a) in newarcs

			#Create a new variable for the arc
			global x[i,a] = @variable(smp, lower_bound = 0) #, upper_bound = 1)
			set_name(x[i,a], string("x[",i,",",a,"]")) 

			#Add to the objective
			set_objective_function(smp, objective_function(smp) + c[a]*x[i,a])

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

		#----------TERMINATION----------#

		if (minimum(minreducedcosts) >= -0.0001) & (setvariables == [])
			println("NO NEGATIVE REDUCED COSTS FOUND!")	
            println("PROCEEDING TO VARIABLE SETTING...")	
			
            varssetfor = [[] for i in orders]
            for i in orders, a in magarcs.A[i]
                if variablefixingthreshold < value(x[i,a]) < 1 - 1e-4 
                    push!(setvariables, (i,a))
                    push!(varssetfor[i], a)
                end
            end
            @constraint(smp, setvars_con[i in orders, a in varssetfor[i]], x[i,a] == 1)

            if variablefixingthreshold == 1.0
                println("NO VARIABLES TO SET")
                break
            end

        elseif (minimum(minreducedcosts) >= -0.0001) & (setvariables != [])
			println("NO NEGATIVE REDUCED COSTS FOUND!")	
			break
		end

		#------------ITERATE------------#

		cg_iter += 1
	
	end

	optimize!(smp)
	smpobj = objective_value(smp)
	#writemagresults(resultsfilename)

    totalarcs = sum(length(magarcs.A[i]) for i in orders)

	#-------------------------------#

	return smpobj, smp, x, y, z, w, magarcs, sum(smptimes), sum(pptimes), sum(pptimes_par), totalarcs

end

