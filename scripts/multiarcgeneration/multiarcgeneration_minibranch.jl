
using Plots
using Base.Iterators: partition

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

function getdualvalues(smpconstraints)

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

### TEST FOR SPEED (MAYBE SCALE UP INSTANCE?)

function findarcvariablereducedcosts(M, smpconstraints)

    alpha, beta, gamma, theta, nu, mu, xi, psi = getdualvalues(smpconstraints)
    
    arcredcosts = zeros(numorders, extendednumarcs)
    for i in orders
        arcredcosts[i,:] += c[1:extendednumarcs] + M.alpha[i] * alpha[i,:] + M.beta[i] * beta[i] + M.gamma[i] * gamma[i] + M.psi[i] * psi[i]
        arcredcosts[i,:] += M.theta * theta + M.nu * nu + M.mu * mu + M.xi * xi
    end

	return arcredcosts

end

#----------------------------------------------------------------------------------------#

function strengthenreducedcosts(arcredcosts, z)

	driverredcost = reduced_cost.(z)

	bestdriver_rc = Dict()
	for a in 1:numarcs
		listoffragments = []
		for d in drivers, l = driverHomeLocs[d], s = drivershift[d], f in fragmentscontaining[l,s,a]
			push!(listoffragments, (d,f))
		end
		if listoffragments != []
			bestdriver_rc[a] = max(0, minimum(driverredcost[frag] for frag in listoffragments))
		else
			bestdriver_rc[a] = 0
		end
	end

	for i in orders, a in 1:numarcs
		arcredcosts[i,a] += bestdriver_rc[a]
	end

	return arcredcosts

end

#----------------------------------------------------------------------------------------#

### TEST FOR SPEED (MAYBE SCALE UP INSTANCE?)

function findarcvariablereducedcosts_varsetting(M, smpconstraints, setvariables)

    alpha, beta, gamma, theta, nu, mu, xi, psi = getdualvalues(smpconstraints)
    
    arcredcosts = zeros(numorders, extendednumarcs)
    for i in orders
        arcredcosts[i,:] += c[1:extendednumarcs] + M.alpha[i] * alpha[i,:] + M.beta[i] * beta[i] + M.gamma[i] * gamma[i] + M.psi[i] * psi[i]
        arcredcosts[i,:] += M.theta * theta + M.nu * nu + M.mu * mu + M.xi * xi
    end

	if varsettingtype == "x"
		for (i,a) in setvariables
        	arcredcosts[i,a] = 0.0 
		end
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

function checkifpathexists(i, myarc, arcset, subproblemsets)

	#Add dummy arcs
	dummyarcs = setdiff(subproblemsets.arclist[i], 1:extendednumarcs)
	arcset = union(arcset, dummyarcs)
	
	#Split into arcs before and after myarcs
	startnode, endnode = subproblemsets.arclookup[i][myarc]
	starttime, endtime = subproblemsets.nodelookup[startnode][2], subproblemsets.nodelookup[endnode][2]
	arcsbefore = sort([a for a in arcset if subproblemsets.nodelookup[subproblemsets.arclookup[i][a][2]][2] <= starttime], by=x->subproblemsets.nodelookup[subproblemsets.arclookup[i][x][1]][2])
	arcsafter = sort([a for a in arcset if subproblemsets.nodelookup[subproblemsets.arclookup[i][a][1]][2] >= endtime], by=x->subproblemsets.nodelookup[subproblemsets.arclookup[i][x][1]][2])

	#Initialize shortest path algorithm 
	currdistance = repeat([999999999.0],outer=[subproblemsets.numnodes])
	currdistance[subproblemsets.dummyorig] = 0

	#Loop over time-space arcs *before* in order of start time
	for a in arcsbefore
		n_end, n_start = subproblemsets.arclookup[i][a][2], subproblemsets.arclookup[i][a][1]
		if currdistance[n_end] > currdistance[n_start] + 0.000001
			currdistance[n_end] = currdistance[n_start]
		end
	end

	#Get distance to myarc
	distancetomyarc = currdistance[startnode]

	#Reset shortest path and ensure we pass through myarc
	currdistance = repeat([999999999.0],outer=[subproblemsets.numnodes])
	currdistance[endnode] = distancetomyarc

	#Loop over time-space arcs *after* in order of start time
	for a in arcsafter  
		n_end, n_start = subproblemsets.arclookup[i][a][2], subproblemsets.arclookup[i][a][1]
		if currdistance[n_end] > currdistance[n_start] + 0.000001
			currdistance[n_end] = currdistance[n_start]
		end
	end

	if currdistance[subproblemsets.dummydest] < 1e-4
		return true
	else
		return false
	end

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

function managecolumnmemory(magarcs, cg_iter, columnmemory, columnmemorylength, variableselected, allremovedarcs)

	removedvars = []

	if cg_iter > columnmemorylength

		#println("--------------------------- COL MGMT ---------------------------")
		
		#Look back a few iterations (defined by columnmemorylength parameter)
		deletefromiter = cg_iter - columnmemorylength

		#Remove arcs added in the chosen iteration that have not been selected since
		for i in orders
			#println("=== ORDER $i ===")
			#println("Arcs added in iteration = $deletefromiter", columnmemory[i,deletefromiter])
			#println("All selected arcs = ", variableselected[i])
			for a in setdiff(columnmemory[i,deletefromiter], variableselected[i])
				if !((i,a) in allremovedarcs)
					remove!(magarcs.A[i], a)
					remove!(magarcs.A_space[i], a)
					n_start, n_end = arcLookup[a]
					remove!(magarcs.A_minus[i, n_end], a)
					remove!(magarcs.A_plus[i, n_start], a)
					push!(removedvars, (i,a))
				end
			end
		end

		#println("---------------------------------------------------------------")
	end

	return magarcs, removedvars

end

#----------------------------------------------------------------------------------------#

function setvarsformag(setvariables, x, z, magarcs)

	if varsettingtype == "x"
		varssetfor = [[] for i in orders]
		varssetcount = 0
		for i in orders, a in magarcs.A[i]
			if variablefixingthreshold[1] + 1e-4 < value(x[i,a]) < variablefixingthreshold[2] - 1e-4 
				push!(setvariables, (i,a))
				push!(varssetfor[i], a)
				varssetcount += 1
			end
		end
		
	elseif varsettingtype == "z"
		varssetfor = [[] for d in drivers]
		varssetcount = 0
		for d in drivers, f in 1:numfragments[driverHomeLocs[d], drivershift[d]]
			if variablefixingthreshold[1] + 1e-4 < value(z[d,f]) < variablefixingthreshold[2] - 1e-4 
				push!(setvariables, (d,f))
				push!(varssetfor[d], f)
				varssetcount += 1
			end
		end
		
	end

	return setvariables, varssetfor, varssetcount

end

#----------------------------------------------------------------------------------------#

function fractionalhistogram(x, z, magarcs, filenamex, filenamez)

	vals = []
	for i in orders, a in magarcs.A[i] 
		if value(x[i,a]) > 1e-4
			push!(vals, value(x[i,a]))
		end
	end
	hist = histogram(vals, bins = 0:0.05:1.1)
	savefig(hist, filenamex)

	vals = []
	for d in drivers, f in 1:numfragments[driverHomeLocs[d], drivershift[d]]
		if value(z[d,f]) > 1e-4
			push!(vals, value(z[d,f]))
		end
	end
	hist = histogram(vals, bins = 0:0.05:1.1)
	savefig(hist, filenamez)

end

#----------------------------------------------------------------------------------------#

function initmagsets(magarcs)
    
    variableusecount = Dict()
	for i in orders, a in magarcs.A[i]
		variableusecount[i,a] = 0
	end
    startercuts = (vars=Dict(), rhs=Dict(), coeff=Dict())
	starterfixedvars = (varsettingiter=[], allvars=[], d=[], f=[], value=[])

    return variableusecount, startercuts, starterfixedvars
end

#----------------------------------------------------------------------------------------#

function combineorderarcsets(arcset, newarcsset)
    
    for i in orders
        newarcs = setdiff(newarcsset.A[i], arcset.A[i])
        for a in newarcs
            push!(arcset.A[i], a)
            n1,n2 = arcLookup[a]
            l1,t1 = nodesLookup[n1]
            l2,t2 = nodesLookup[n2]
            if !(l1 == l2)
                push!(arcset.A_space[i], a)
            end
            push!(arcset.A_plus[i,n1], a)
            push!(arcset.A_minus[i,n2], a)
        end
    end

    return arcset

end

#----------------------------------------------------------------------------------------#

function multiarcgeneration_minibranch!(magarcs, hasdriverarcs, startercuts, starterfixedvars, variableusecount, currvarfixingiter, cg_iter_start)

	#currvarfixingiter, cg_iter_start = 0, 1

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

    if symmetrybreaking_flag == 1
		for l in 1:numlocs, s in 1:numeffshifts
			for j in 1:length(driversets[l,s])-1
				d1, d2 = driversets[l,s][j], driversets[l,s][j+1]
				@constraint(ip, sum(z[d1,f] for f in setdiff(1:numfragments[l,s], workingfragments[l,s])) <= sum(z[d2,f] for f in setdiff(1:numfragments[l,s], workingfragments[l,s])) )
			end
		end
	end

    #Add any existing cuts
    @constraint(smp, [i in 1:length(startercuts.vars)], sum(startercuts.coeff[i][d,j] * z[d,j] for (d,j) in startercuts.vars[i]) <= startercuts.rhs[i])

    #Fix any existing variables     
    @constraint(smp, [i in starterfixedvars.allvars], z[starterfixedvars.d[i], starterfixedvars.f[i]] == starterfixedvars.value[i])

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
	cg_iter = cg_iter_start
	smpobjectives, smptimes, pptimes, pptimes_par, lowerbounds = [], [], [], [], []
	listlength = convert(Int64, ceil(length(orders)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)

	#Master list of knapsack cuts and fixed variables
	mastercuts = startercuts #(vars=[], coeff=[], rhs=[])

	#------------------------------------------------------#

	#Pre-processing
	println("Starting...")
    M, subproblemsets = preprocessmagsets(orderarcs.A);
	println("Length of subproblemsets.nodelookup = ", length(subproblemsets.nodelookup))

	#------------------------------------------------------#

    setvariables = []
    variableselected = Dict()
    for i in orders
        variableselected[i] = []
    end
	columnmemory, allremovedarcs = Dict(), []
	addcutsthisiter_flag = 0
	cuttime = 0
    bestlowerbound = 0
	arcsbeforevarsetting = (A=Dict(), A_space=Dict())
	for i in orders
		arcsbeforevarsetting.A[i] = []
	end

    #------------------------------------------------------#

	while cg_iter <= 100000

		#-----------SOLVE SMP-----------#

		println("-------- ITERATION $cg_iter --------")             
		status = optimize!(smp)
		if termination_status(smp) != MOI.OPTIMAL
			println(termination_status(smp))
			return 100000000, smp, x, y, z, w, magarcs, 0, 0, 0, 0, 0, (vars=Dict(), rhs=Dict(), coeff=Dict()), 0
		end
		smpobj, smptime = objective_value(smp), solve_time(smp)
    
        #Update chosen variables
        for i in orders, a in magarcs.A[i]
            if value(x[i,a]) > 1e-4
                push!(variableselected[i], a)
				variableusecount[i,a] += max(1,cg_iter/20) #More weight to later iterations
            end
        end
		
		push!(smpobjectives, copy(smpobj))
		push!(smptimes, copy(smptime))	
		println("Solved SMP with objective = ", smpobj, " in iteration $cg_iter (", sum(length(magarcs.A[i]) for i in orders), " arcs)")

        #------------SUBPROBLEMS------------#

		#Calculate reduced costs
        if setvariables == []
            arcredcosts = findarcvariablereducedcosts(M, smpconstraints)
        else
            arcredcosts = findarcvariablereducedcosts_varsetting(M, smpconstraints, setvariables)
        end
		if strengthenedreducedcost_flag == 1
			arcredcosts = strengthenreducedcosts(arcredcosts, z)
		end

		#Run shortest path for each order to find new arcs
		dptimelist = []
		addarcs, minreducedcosts = [], []
		for i in orders
			if solutionmethod == "mag"	
				minreducedcost, shortestpathnodes, shortestpatharcs, sptime = magsubproblem(i, arcredcosts, subproblemsets)
			elseif solutionmethod == "sag"	
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

		#-----------ADD CUTS-----------#
	
		if addcutsthisiter_flag == 1
			cutstarttime = time()
			cuts = findknapsackcuts(z, knapsackcuttype)
			@constraint(smp, [i in 1:cuts.numcuts], sum(cuts.coeff[i][d,j] * z[d,j] for (d,j) in cuts.vars[i]) <= cuts.rhs[i])
			cutindex = length(mastercuts.rhs)+1
			for i in 1:cuts.numcuts
				mastercuts.vars[cutindex] = cuts.vars[i]
				mastercuts.coeff[cutindex] = cuts.coeff[i]
				mastercuts.rhs[cutindex] = cuts.rhs[i]
				cutindex += 1
			end
			cutsaddedthisiter = cuts.numcuts
			cuttime += time() - cutstarttime
			println("Added ", cuts.numcuts, " cuts")
		else 
			cutsaddedthisiter = -1 * knapsackcuts_flag
		end

		#-------ADD NEW VARIABLES-------#

		#Add new arcs to order arc sets
		magarcs, newarcs = updatearcsets(magarcs, addarcs)

		#Add new row to the column memory
		for i in orders
			columnmemory[i,cg_iter] = []
		end

		#Add new arcs to model as x-variables
		for (i,a) in newarcs

			push!(columnmemory[i,cg_iter], a)
			variableusecount[i,a] = currvarfixingiter*100 #Never remove arcs added via variable fixing

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

        if (currvarfixingiter == 0) & (minimum(minreducedcosts) >= -0.0001) & (cutsaddedthisiter == 0)
            bestlowerbound += smpobj
        end

		if (minimum(minreducedcosts) >= -0.0001) & (cutsaddedthisiter == -1) 
			addcutsthisiter_flag += 1
		elseif (minimum(minreducedcosts) >= -0.0001) & (currvarfixingiter < varsettingiterations) & (cutsaddedthisiter == 0)
			println("NO NEGATIVE REDUCED COSTS FOUND!")	
            println("PROCEEDING TO VARIABLE SETTING...")

			if currvarfixingiter == 0
				for i in orders, a in magarcs.A[i]
					push!(arcsbeforevarsetting.A[i], a)
				end 
			end
			
			#fractionalhistogram(x, z, magarcs, string("outputs/activevars_x_", currvarfixingiter,".png"), string("outputs/activevars_z_", currvarfixingiter,".png"))
           
            #Select the variables to fix based on pre-defined thresholds
            setvariables, varssetfor, varssetcount = setvarsformag(setvariables, x, z, magarcs)
            
            #If no variables are within the thresholds, then this round of MAG is complete 
            if (variablefixingthreshold[2] == 1.0) || (varssetcount == 0)
                println("NO VARIABLES TO SET")
                break
            end

            #Kick-off two new MAG nodes, one setting the fixed vars to 1 and the other to 0
            for fixedvalue in [1,0]

                #Format the fixed variable info to pass along to the new node
                masterfixedvars = deepcopy(starterfixedvars) 
                for d in drivers, f in varssetfor[d]
                    push!(masterfixedvars.allvars, length(masterfixedvars.allvars)+1)
                    push!(masterfixedvars.d, d)
                    push!(masterfixedvars.f, f)
                    push!(masterfixedvars.value, fixedvalue)
                end
                push!(masterfixedvars.varsettingiter, fixedvalue)

                println("Variables used before var setting = ", sum(values(variableusecount)))

                #Run MAG with the new fixed variables
                println("------------------- SETTING VARIABLES - ", masterfixedvars.varsettingiter, " -------------------")
                obj_fix, smp_fix, x_fix, y_fix, z_fix, w_fix, magarcs_fix, smptm_fix, pptm_fix, pppar_fix, arcs_fix, cgiter_fix, cuts_fix, cuttime = multiarcgeneration_minibranch!(magarcs, hasdriverarcs, mastercuts, masterfixedvars, variableusecount, currvarfixingiter+1, cg_iter+1)
                cg_iter += - cg_iter + cgiter_fix 

                println("Variables used after var setting = ", sum(values(variableusecount)))

                #Add the generated arcs to the master list
                println("Arcs before var setting = ", sum(length(magarcs.A[i]) for i in orders))
                magarcs = combineorderarcsets(magarcs, magarcs_fix)
                println("Arcs after var setting = ", sum(length(magarcs.A[i]) for i in orders))
                
                #Add the generated cuts to the master list
                println("Cuts before var setting = ", length(mastercuts.rhs))
                cutindex = length(mastercuts.rhs)+1
                for i in 1:length(cuts_fix.rhs)
                    mastercuts.vars[cutindex] = cuts_fix.vars[i]
                    mastercuts.coeff[cutindex] = cuts_fix.coeff[i]
                    mastercuts.rhs[cutindex] = cuts_fix.rhs[i]
                    cutindex += 1
                end
                println("Cuts after var setting = ", length(mastercuts.rhs))
                                    
            end

            println("COMPLETED ALL VAR SETTING NODES")
            break

        elseif (minimum(minreducedcosts) >= -0.0001) & (currvarfixingiter >= varsettingiterations) & (cutsaddedthisiter == 0)
			println("NO NEGATIVE REDUCED COSTS FOUND!")	
			break
		elseif (minimum(minreducedcosts) >= -0.0001) & (cutsaddedthisiter >= 1)
			println("Keep adding cuts")
		end

		#------------ITERATE------------#

		cg_iter += 1
	
	end

	arcsbeforecolmgmt = (A=[copy(magarcs.A[i]) for i in orders], A_space=Dict())

    if currvarfixingiter == 0
        
        println("Before arcs = ", sum(length(magarcs.A[i]) for i in orders))
        
        #Column management
        for i in orders, a in intersect(1:numarcs, magarcs.A[i])
            if variableusecount[i,a] / 100 <= postmagcolumndeletionthreshold
                remove!(magarcs.A[i], a)
                remove!(magarcs.A_space[i], a)
                n_start, n_end = arcLookup[a]
                remove!(magarcs.A_minus[i, n_end], a)
                remove!(magarcs.A_plus[i, n_start], a)
            end
        end

		println("After arcs = ", sum(length(magarcs.A[i]) for i in orders))

		#Clear out any unusable arcs
		for i in orders, a in intersect(1:numarcs, magarcs.A[i])
			if !(checkifpathexists(i, a, magarcs.A[i], subproblemsets))
				remove!(magarcs.A[i], a)
                remove!(magarcs.A_space[i], a)
                n_start, n_end = arcLookup[a]
                remove!(magarcs.A_minus[i, n_end], a)
                remove!(magarcs.A_plus[i, n_start], a)
			end
		end

        println("After after arcs = ", sum(length(magarcs.A[i]) for i in orders))

    end
        
    #-------------------------------#

    totalarcs = sum(length(magarcs.A[i]) for i in orders)

	#-------------------------------#

    return bestlowerbound, smp, x, y, z, w, magarcs, sum(smptimes), sum(pptimes), sum(pptimes_par), totalarcs, cg_iter, mastercuts, cuttime, arcsbeforecolmgmt, arcsbeforevarsetting
    
end


