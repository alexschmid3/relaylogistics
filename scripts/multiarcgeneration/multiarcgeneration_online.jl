
using Plots
using Base.Iterators: partition

include("preprocessmagsets_online.jl")

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

    alpha = converttosparsearray(dual.(smpconstraints.con_orderFlowBalance), maximum(currstate.orders), extendednumnodes)
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
	for i in currstate.orders, a in 1:extendednumarcs
		arcredcosts[i,a] = c[a] 
	end
    
	#Calculate reduced costs for each arc
	for i in currstate.orders, n in setdiff([n2 for n2 in 1:numnodes], union(currstate.Origin[i], currstate.Destination[i])), a in A_minus[n]
		arcredcosts[i,a] -= alpha[i,n]
	end
	for i in currstate.orders, n in setdiff([n2 for n2 in 1:numnodes], union(currstate.Origin[i], currstate.Destination[i])), a in A_plus[n]
		arcredcosts[i,a] += alpha[i,n] 
	end

	for i in currstate.orders, n in currstate.Origin[i], a in intersect(union(A_space, dummyarc), A_plus[n])
		arcredcosts[i,a] -= beta[i] 
	end
	for i in setdiff(currstate.orders, currstate.ordersinprogress)
		a = extendedarcs[last(currstate.Origin[i]), last(currstate.Destination[i])]
		arcredcosts[i,extendedarcs[last(currstate.Origin[i]), last(currstate.Destination[i])]] -= beta[i] 
	end
	
	for i in intersect(currstate.orders, currstate.ordersinprogress), n in currstate.Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
		arcredcosts[i,a] -= beta[i] 
	end
	
	for i in currstate.orders, n in currstate.Destination[i], a in A_minus[n]
		if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))
			arcredcosts[i,a] -= gamma[i]
		end
	end
	
	for n in N_0, i in currstate.orders, a in setdiff(A_plus[n], dummyarc)  
		arcredcosts[i,a] -= theta[n]
	end
	for n in N_end, a in A_minus[n], i in currstate.orders
		arcredcosts[i,a] -= nu[n]
	end
	for n in N_flow_t, a in A_minus[n], i in currstate.orders
		arcredcosts[i,a] -= mu[n]
	end
	for n in N_flow_t, a in A_plus[n], i in currstate.orders
		arcredcosts[i,a] += mu[n]
	end
	
	for i in currstate.orders, n in currstate.Destination[i], a in A_minus[n]
		arcredcosts[i,a] += arcfinishtime[a] * psi[i]
	end
	
    for i in currstate.orders, a in A_space
		arcredcosts[i,a] += xi[a]
	end
    

	return arcredcosts

end

#----------------------------------------------------------------------------------------#

### TEST FOR SPEED (MAYBE SCALE UP INSTANCE?)

function findarcvariablereducedcosts(M, smpconstraints)

    alpha, beta, gamma, theta, nu, mu, xi, psi = getdualvalues(smpconstraints)
    
    arcredcosts = zeros(maximum(currstate.orders), extendednumarcs)
    for i in currstate.orders
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

	for i in currstate.orders, a in 1:numarcs
		arcredcosts[i,a] += bestdriver_rc[a]
	end

	return arcredcosts

end

#----------------------------------------------------------------------------------------#

function sagsubproblem(i, arcredcosts)

    spstarttime = time()

    mostnegativearc = setdiff(currarcs.orderarcs.A[i], dummyarc)[argmin([arcredcosts[i,a] for a in setdiff(currarcs.orderarcs.A[i], dummyarc)])]
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

function updatearcsets(currarcs, addarcs)

    newarcs = []

    for (i,a) in addarcs
        if !(a in currarcs.magarcs.A[i])
            push!(currarcs.magarcs.A[i], a)
            if a in primaryarcs.A_space
                push!(currarcs.magarcs.A_space[i], a)
            end
            n_plus = arcLookup[a][1]
            n_minus = arcLookup[a][2]
            push!(currarcs.magarcs.A_plus[i,n_plus], a)
            push!(currarcs.magarcs.A_minus[i,n_minus], a)
            push!(newarcs, (i,a))
        end
    end

    return currarcs, newarcs

end

#----------------------------------------------------------------------------------------#

function multiarcgeneration!(currstate, currfragments, currarcs)

    smp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(smp, "TimeLimit", 60*60)
	set_optimizer_attribute(smp, "OutputFlag", 0)
    set_optimizer_attribute(smp, "MIPGap", 0.001)

	#Variables
	x = Dict()
	for i in currstate.orders, a in currarcs.magarcs.A[i]
	    global x[i,a] = @variable(smp, lower_bound = 0) #, upper_bound = 1) <-- UB implied by order origin/dest con & messes up reduced cost sign
	    set_name(x[i,a], string("x[",i,",",a,"]")) 
	end
	@variable(smp, y[currarcs.hasdriverarcs.A] >= 0)
    @variable(smp, 0 <= z[d = drivers, f = 1:currfragments.numfragments[driverHomeLocs[d], drivershift[d]]] <= 1)
	@variable(smp, w[a in primaryarcs.A_space] >= 0)
	@variable(smp, ordtime[currstate.orders])

	#Objective
	@objective(smp, Min, lambda * sum((ordtime[i] - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i] for i in currstate.orders) + sum(sum(c[a]*x[i,a] for a in currarcs.magarcs.A[i]) for i in currstate.orders) + sum(c[a]*y[a] for a in currarcs.hasdriverarcs.A) + sum(u[a]*w[a] for a in primaryarcs.A_space))
		
	#Order constraints
	@constraint(smp, orderFlowBalance[i = currstate.orders, n in setdiff([n2 for n2 in 1:numnodes], union(currstate.Origin[i], currstate.Destination[i]))], sum(x[i,a] for a in currarcs.magarcs.A_minus[i,n]) - sum(x[i,a] for a in currarcs.magarcs.A_plus[i,n]) == 0)
	@constraint(smp, arriveDestin[i = currstate.orders], sum(sum(x[i,a] for a in currarcs.magarcs.A_minus[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in currstate.Destination[i]) == 1)
	@constraint(smp, departOrigin[i = currstate.orders], sum(sum(x[i,a] for a in intersect(union(primaryarcs.A_space, dummyarc), currarcs.magarcs.A_plus[i,n])) for n in currstate.Origin[i]) == 1)
	for i in setdiff(currstate.orders, currstate.ordersinprogress)
		extendedorderarc = extendedarcs[last(currstate.Origin[i]), last(currstate.Destination[i])]
		if extendedorderarc in currarcs.magarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end
	
	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(currstate.orders, currstate.ordersinprogress), n in currstate.Origin[i], a in setdiff(primaryarcs.A_plus[n], union(primaryarcs.A_space, dummyarc))
		if a in currarcs.magarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(smp, deliveryTime[i in currstate.orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in currarcs.magarcs.A_minus[i,n]) for n in currstate.Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(smp, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(currarcs.magarcs.A_plus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == currstate.m_0[n])
	@constraint(smp, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(currarcs.magarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) >= currstate.m_end[n])
	@constraint(smp, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(currarcs.magarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n])- sum(sum(x[i,a] for a in setdiff(currarcs.magarcs.A_plus[i,n],dummyarc)) for i in currstate.orders) - sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == 0)

	#Linking constraints
	@constraint(smp, driverAvailability[a in primaryarcs.A_space], sum(sum(z[d,f] for f in currfragments.fragmentscontaining[driverHomeLocs[d], drivershift[d],a]) for d in drivers) == w[a]  )
	for i in currstate.orders, a in currarcs.magarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in currarcs.hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
    @constraint(smp, driverStartingLocs[d in drivers, l = driverHomeLocs[d], s = drivershift[d]], sum(sum(z[d,f] for f in currfragments.F_plus_ls[l,s,n]) for n in currfragments.driverSetStartNodes[l,s]) == 1)
	@constraint(smp, driverFlowBalance[d in drivers, l = driverHomeLocs[d], s = drivershift[d], n in currfragments.N_flow_ls[l,s]], sum(z[d,f] for f in currfragments.F_minus_ls[l,s,n]) - sum(z[d,f] for f in currfragments.F_plus_ls[l,s,n]) == 0)

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
	listlength = convert(Int64, ceil(length(currstate.orders)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)

	#------------------------------------------------------#

	#Pre-processing
	println("Starting...")
    M, subproblemsets = preprocessmagsets(currarcs.orderarcs.A, currstate);
	println("Length of subproblemsets.nodelookup = ", length(subproblemsets.nodelookup))

	#------------------------------------------------------#

    variableselected = Dict()
    for i in currstate.orders
        variableselected[i] = []
    end
	columnmemory, allremovedarcs = Dict(), []
    bestlowerbound = 0

    #------------------------------------------------------#

	while cg_iter <= 100000

		#-----------SOLVE SMP-----------#

		println("-------- ITERATION $cg_iter --------")             
		status = optimize!(smp)
		if termination_status(smp) != MOI.OPTIMAL
			println(termination_status(smp))
			return 100000000, smp, x, y, z, w, currarcs.magarcs, 0, 0, 0, 0, 0, (vars=Dict(), rhs=Dict(), coeff=Dict()), 0
		end
		smpobj, smptime = objective_value(smp), solve_time(smp)
		push!(smpobjectives, copy(smpobj))
		push!(smptimes, copy(smptime))	
		println("Solved SMP with objective = ", smpobj, " in iteration $cg_iter (", sum(length(currarcs.magarcs.A[i]) for i in currstate.orders), " arcs)")

        #------------SUBPROBLEMS------------#

		#Calculate reduced costs
        arcredcosts = findarcvariablereducedcosts(M, smpconstraints)
		if strengthenedreducedcost_flag == 1
			arcredcosts = strengthenreducedcosts(arcredcosts, z)
		end

		#Run shortest path for each order to find new arcs
		dptimelist = []
		addarcs, minreducedcosts = [], []
		for i in currstate.orders
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
		shuffleddptimes = shuffle_partition(length(currstate.orders))
		dptimelistsums = [sum(dptimelist[j_inner] for j_inner in shuffleddptimes[j_outer]) for j_outer in 1:length(shuffleddptimes)]
		push!(pptimes_par, maximum(dptimelistsums))
        push!(pptimes, sum(dptimelist))
		
        #Update the lowerbound 
        try
			push!(lowerbounds, smpobj + sum(minreducedcosts[k] for k in 1:length(minreducedcosts) if minreducedcosts[k] <= -0.000001))
		catch
			push!(lowerbounds, smpobj)
		end

		#-------ADD NEW VARIABLES-------#

		#Add new arcs to order arc sets
		currarcs, newarcs = updatearcsets(currarcs, addarcs)

		#Add new row to the column memory
		for i in currstate.orders
			columnmemory[i,cg_iter] = []
		end

		#Add new arcs to model as x-variables
		for (i,a) in newarcs

			#Create a new variable for the arc
			global x[i,a] = @variable(smp, lower_bound = 0) #, upper_bound = 1)
			set_name(x[i,a], string("x[",i,",",a,"]")) 

			#Add to the objective
			set_objective_coefficient(smp, x[i,a], c[a])

			#Find the entering and exiting nodes for arc a (ex. n_plus is the node for which a belongs to A_plus[n])
			n_plus = arcLookup[a][1]
			n_minus = arcLookup[a][2]

			#Add new variable to order constraints
			if n_minus in setdiff([n2 for n2 in 1:numnodes], union(currstate.Origin[i], currstate.Destination[i]))
				set_normalized_coefficient(smpconstraints.con_orderFlowBalance[i,n_minus], x[i,a], 1.0)
			end
			if n_plus in setdiff([n2 for n2 in 1:numnodes], union(currstate.Origin[i], currstate.Destination[i]))
				set_normalized_coefficient(smpconstraints.con_orderFlowBalance[i,n_plus], x[i,a], -1.0)
			end
			if (n_plus in currstate.Origin[i]) & (a in A_space)
				set_normalized_coefficient(smpconstraints.con_departOrigin[i], x[i,a], 1.0)
			end
			if (n_plus in currstate.Origin[i]) & (i in intersect(currstate.orders, currstate.ordersinprogress)) & !(a in union(A_space, dummyarc))
				set_normalized_coefficient(smpconstraints.con_departOrigin[i], x[i,a], 1.0)
			end
			if (inprogressdummyarc_flag == 0) && !(i in currstate.ordersinprogress) && (a == extendedarcs[last(currstate.Origin[i]), last(currstate.Destination[i])]) 
				set_normalized_coefficient(smpconstraints.con_departOrigin[i], x[i,a], 1.0)
			end
			if (n_minus in currstate.Destination[i]) & ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))
				set_normalized_coefficient(smpconstraints.con_arriveDestin[i], x[i,a], 1.0)
			end

			#Add new variable to order delivery constraints
			if (n_minus in currstate.Destination[i]) & (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1])
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

		if (minimum(minreducedcosts) >= -0.0001) 
			println("NO NEGATIVE REDUCED COSTS FOUND!")	
			break
		end

		#------------ITERATE------------#

		cg_iter += 1
	
	end
        
    #-------------------------------#

    totalarcs = sum(length(currarcs.magarcs.A[i]) for i in currstate.orders)

	#-------------------------------#

    return bestlowerbound, smp, x, y, z, w, currarcs.magarcs, sum(smptimes), sum(pptimes), sum(pptimes_par), totalarcs
    
end


