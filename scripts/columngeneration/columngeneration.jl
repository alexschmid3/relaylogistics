
include("preprocesscgsets.jl")

function initialfeasiblepaths(orderarcs)
		
	#--------------------------CREATE --------------------------#

	#Create delta[i,a,p] = 1 if path p for order i contains arc a, 0 otherwise
	delta = Dict()

	#Create paths[i] = list of path indices for order i (indices added with column generation)
	paths = Dict()
	for i in orders
		paths[i] = [1, 2]
	end

	#--------------------------ADD DUMMY PATHS--------------------------#

	#Initialize delta
	for i in orders, a in orderarcs.A[i], p in paths[i]
		delta[i,a,p] = 0
	end
	for i in orders
		delta[i, dummyarc, 1] = 1
		delta[i, extendedarcs[nodes[originloc[i], horizon], extendednodes[destloc[i], dummyendtime]], 2] = 1
	end
	
	return delta, paths

end

#----------------------------------------------------------------------------------------#

function getdualvalues_cg(rmpconstraints)

	alpha = dual.(rmpconstraints.con_orderpath)
	beta = Array(dual.(rmpconstraints.con_initialTrucks))
	epsilon = Array(dual.(rmpconstraints.con_finalTrucks))
	gamma = Array(dual.(rmpconstraints.con_truckFlowBalance))
	eta = dual.(rmpconstraints.con_driverAvailability)
	psi = dual.(rmpconstraints.con_deliveryTime)

    return alpha, beta, epsilon, gamma, eta, psi

end

#----------------------------------------------------------------------------------------#

### TEST FOR SPEED (MAYBE SCALE UP INSTANCE?)

function findarcvariablereducedcosts_cg(M, rmpconstraints)

    alpha, beta, epsilon, gamma, eta, psi = getdualvalues_cg(rmpconstraints)
    
    arcredcosts = zeros(numorders, extendednumarcs)
    for i in orders
        arcredcosts[i,:] += c[1:extendednumarcs] + M.psi[i] * psi[i]
        arcredcosts[i,:] += M.beta * beta + M.eta * eta + M.epsilon * epsilon + M.gamma * gamma
    end

	return arcredcosts, alpha

end

#----------------------------------------------------------------------------------------#

function cgsubproblem(i, arcredcosts, subproblemsets)

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

function columngeneration!(orderarcs, hasdriverarcs)
	
	#Initialize paths
	delta, paths = initialfeasiblepaths(orderarcs)

    #------------------------------------------------------#
	
	# Restricted master problem
	rmp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(rmp, "TimeLimit", 60*60)
	set_optimizer_attribute(rmp, "OutputFlag", 0)
    set_optimizer_attribute(rmp, "MIPGap", 0.001)

	#Variables
	x = Dict()
	for i in orders, p in paths[i] 
	    global x[i,p] = @variable(rmp, lower_bound = 0) #, upper_bound = 1)
	    set_name(x[i,p], string("x[",i,",",p,"]")) 
	end
	@variable(rmp, y[hasdriverarcs.A] >= 0)
    @variable(smp, 0 <= z[d = drivers, f = 1:numfragments[driverHomeLocs[d], drivershift[d]]] <= 1)
	@variable(rmp, w[a in primaryarcs.A_space] >= 0)
	@variable(rmp, ordtime[orders])
    @variable(rmp, maxhours)

	#Objective
	@objective(rmp, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(sum(c[a]*delta[i,a,p]*x[i,p] for a in orderarcs.A[i] ) for p in paths[i]) for i in orders) + sum(c[a]*y[a] for a in hasdriverarcs.A) + sum(u[a]*w[a] for a in primaryarcs.A_space) + lambda2 * maxhours)
	
	#Order constraints
	@constraint(rmp, orderpath[i in orders], sum(x[i,p] for p in paths[i]) == 1)
	
	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(orders, ordersinprogress), n in Origin[i], a in setdiff(primaryarcs.A_plus[n], union(primaryarcs.A_space, dummyarc))
		if a in orderarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(rmp, deliveryTime[i in orders], ordtime[i] - sum(sum(sum(arcfinishtime[a] * delta[i,a,p] * x[i,p] for a in orderarcs.A_minus[i,n]) for p in paths[i]) for n in Destination[i]) == - orderOriginalStartTime[i] )

	#Truck constraints
	@constraint(rmp, initialTrucks[n in N_0], sum(sum(sum(delta[i,a,p] * x[i,p] for a in orderarcs.A_plus[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_plus[n]) == m_0[n])
	@constraint(rmp, finalTrucks[n in N_end], sum(sum(sum(delta[i,a,p] * x[i,p] for a in orderarcs.A_minus[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) >= m_end[n])
	@constraint(rmp, truckFlowBalance[n in N_flow_t], sum(sum(sum(delta[i,a,p] * x[i,p] for a in orderarcs.A_minus[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) - sum(sum(sum(delta[i,a,p] * x[i,p] for a in orderarcs.A_plus[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) - sum(y[a] for a in hasdriverarcs.A_plus[n]) == 0)

	#Linking constraints
	@constraint(rmp, driverAvailability[a in primaryarcs.A_space], sum(sum(z[d,f] for f in fragmentscontaining[driverHomeLocs[d], drivershift[d],a]) for d in drivers) == w[a] )
	for i in orders, p in paths[i], a in orderarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,p], -1 * delta[i,a,p])
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
	rmpconstraints = (con_orderpath = orderpath,
		con_initialTrucks = initialTrucks,
		con_finalTrucks = finalTrucks,
		con_truckFlowBalance = truckFlowBalance,
		con_driverAvailability = driverAvailability,
		con_deliveryTime = deliveryTime
		)

	#------------------------------------------------------#

	#Initialize column generation 
	cg_iter = 1
	rmpobjectives, rmptimes, pptimes, pptimes_par, lowerbounds = [], [], [], [], []
	listlength = convert(Int64, ceil(length(orders)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)

	#------------------------------------------------------#

	#Pre-processing
    M, subproblemsets = preprocesscgsets(orderarcs.A);

    #------------------------------------------------------#

	while cg_iter <= 100000

		#-------------SOLVE RMP-------------#

		println("-------- ITERATION $cg_iter --------")
		status = optimize!(rmp)
		if termination_status(rmp) != MOI.OPTIMAL
			println(termination_status(rmp))
			return 100000000, rmp, x, y, z, w, paths, delta
		end
		rmpobj, rmptime = objective_value(rmp), solve_time(rmp)
		push!(rmpobjectives, copy(rmpobj))
		push!(rmptimes, copy(rmptime))	
		println("Solved RMP with objective = ", rmpobj, " in iteration $cg_iter (", sum(length(paths[i]) for i in orders), " paths)")

		#Count number of arcs and paths
		if saveconvergencedata_flag >= 0
			totalorderarcs = sum(sum(min(1, sum(delta[i,a,p] for p in paths[i])) for a in orderarcs.A[i]) for i in orders)
			totalorderpaths = sum(length(paths[i]) for i in orders)
			totalusedpaths = sum(sum(value(x[i,p]) for p in paths[i]) for i in orders)
		end

        #------------SUBPROBLEMS------------#

		#Calculate reduced costs
        arcredcosts, alpha = findarcvariablereducedcosts_cg(M, rmpconstraints)

		#Run shortest path for each order to find new arcs
		shortestpathnodes, shortestpatharcs = Dict(), Dict()
		dptimelist, newpaths, minreducedcosts,  = [], [], []
		for i in orders
			pathredcost, spnodes, sparcs, sptime = cgsubproblem(i, arcredcosts, subproblemsets)
			
			#Add the dual value for the "one order, one path" constraint
			minreducedcost = pathredcost - alpha[i]

			push!(dptimelist, sptime)
			push!(minreducedcosts, minreducedcost)

			if minreducedcost < -0.0001
				newpathindex = last(paths[i]) + 1
				pathobjcost = 0
				for a in orderarcs.A[i]
					delta[i,a,newpathindex] = 0
				end
				for a in sparcs
					delta[i,a,newpathindex] = 1
					pathobjcost += c[a]
				end
				push!(paths[i], newpathindex)
				push!(newpaths, (i, newpathindex, pathobjcost))
				shortestpathnodes[i] = spnodes
				shortestpatharcs[i] = sparcs			
			end
		end
		
		#"Parallelize" subproblem times
		shuffleddptimes = shuffle_partition(length(orders))
		dptimelistsums = [sum(dptimelist[j_inner] for j_inner in shuffleddptimes[j_outer]) for j_outer in 1:length(shuffleddptimes)]
		push!(pptimes_par, maximum(dptimelistsums))
        push!(pptimes, sum(dptimelist))
		
        #Update the lowerbound 
        try
			push!(lowerbounds, rmpobj + sum(minreducedcosts[k] for k in 1:length(minreducedcosts) if minreducedcosts[k] <= -0.000001))
		catch
			push!(lowerbounds, rmpobj)
		end

		#------COUNT ARCS AND PATHS-----#

		if saveconvergencedata_flag >= 0
			maximprove = minimum(minreducedcosts) * totalusedpaths
			write_cg_conv(convergencedatafilename, cg_iter, maximprove, totalorderarcs, totalorderpaths, rmpobj)
		end
	
		#-------ADD NEW VARIABLES-------#

		for row in newpaths
			#Get order index, path index, and path cost
			i, p = row[1], row[2]
			cost = sum(c[a]*delta[i,a,p] for a in orderarcs.A[i]) #row[3]
		
			#Create a new variable for the path
			global x[i,p] = @variable(rmp, upper_bound = 1, lower_bound = 0)
			set_name(x[i,p], string("x[",i,",",p,"]")) 

			#Add to objective
			set_objective_function(rmp, objective_function(rmp) + cost*x[i,p]) 

			#Add new variable to order path constraints
			set_normalized_coefficient(rmpconstraints.con_orderpath[i], x[i,p], 1.0)

			#Add to delivery time constraintsinclude("pathbasedcolumngeneration.jl")
			set_normalized_coefficient(rmpconstraints.con_deliveryTime[i], x[i,p], - sum(sum(arcfinishtime[a] * delta[i,a,p] for a in orderarcs.A_minus[i,n]) for n in Destination[i] if orderarcs.A_minus[i,n] != []) )
			
			#Add new variable to initial and final trucks constraints
			for l in 1:numlocs
				if orderarcs.A_plus[i,nodes[l,0]] != []
					set_normalized_coefficient(rmpconstraints.con_initialTrucks[nodes[l,0]], x[i,p], sum(delta[i,a,p] for a in orderarcs.A_plus[i,nodes[l,0]]))
				end
				if orderarcs.A_minus[i,nodes[l,horizon]] != []
					set_normalized_coefficient(rmpconstraints.con_finalTrucks[nodes[l,horizon]], x[i,p], sum(delta[i,a,p] for a in orderarcs.A_minus[i,nodes[l,horizon]]))
				end
			end

			#Add new variable to truck flow balance constraints
			n = minimum(shortestpathnodes[i])
			if n in N_flow_t
				set_normalized_coefficient(rmpconstraints.con_truckFlowBalance[n], x[i,p], -1.0 )
			end
			n = maximum(shortestpathnodes[i])
			if n in N_flow_t
				set_normalized_coefficient(rmpconstraints.con_truckFlowBalance[n], x[i,p], 1.0 )
			end

			#Add new variable to driver availability constraints
			for a in intersect(orderarcs.A_space[i], shortestpatharcs[i])
				set_normalized_coefficient(rmpconstraints.con_driverAvailability[a], x[i,p], -1.0)
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

	optimize!(rmp)
	rmpobj = objective_value(rmp)

    #-------------------------------#

	totalpaths = sum(length(paths[i]) for i in orders)

	#-------------------------------#

	return rmpobj, rmp, x, y, z, w, paths, delta, sum(rmptimes), sum(pptimes), sum(pptimes_par), totalpaths, cg_iter

end

