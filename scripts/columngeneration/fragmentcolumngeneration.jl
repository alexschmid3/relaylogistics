
include("convergenceplots.jl")

using Base.Iterators: partition

#----------------------------------------------------------------------------------------#

function initialfeasiblesolution_pbcg(orders, orderArcSet)
		
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
	for i in orders, a in orderArcSet[i], p in paths[i]
		delta[i,a,p] = 0
	end
	for i in orders
		delta[i, dummyarc, 1] = 1
		delta[i, extendedarcs[nodes[originloc[i], horizon], extendednodes[destloc[i], dummyendtime]], 2] = 1
	end
	
	return delta, paths

end

#----------------------------------------------------------------------------------------#

function sparsemasterproblem_pbcg(paths, delta, orderArcSet_full, orderArcSet_space_full, numeffshifts)

	rmp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(rmp, "TimeLimit", 60*60)
	set_optimizer_attribute(rmp, "OutputFlag", 0)

	#Variables
	x = Dict()
	for i in orders, p in paths[i] 
	    global x[i,p] = @variable(rmp, lower_bound = 0) #, upper_bound = 1)
	    set_name(x[i,p], string("x[",i,",",p,"]")) 
	end
	@variable(rmp, y[A_hasdriver] >= 0)
	@variable(rmp, z[l = 1:numlocs, s = 1:numeffshifts, f = 1:numfragments[l,s]] >= 0)
	#@variable(rmp, 0 <= z[d = drivers, homeArcSet[d]] <= 1)
	@variable(rmp, w[a in A_space] >= 0)
	@variable(rmp, ordtime[orders])

	#Objective
	@objective(rmp, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(sum(c[a]*delta[i,a,p]*x[i,p] for a in orderArcSet_full[i] ) for p in paths[i]) for i in orders) + sum(c[a]*y[a] for a in A_hasdriver) + sum(u[a]*w[a] for a in A_space) )

	#Order constraints
	@constraint(rmp, orderpath[i in orders], sum(x[i,p] for p in paths[i]) == 1)

	#Order delivery constraints
	@constraint(rmp, deliveryTime[i in orders], ordtime[i] - sum(sum(sum(arcfinishtime[a] * delta[i,a,p] * x[i,p] for a in A_minus_i_full[i,n]) for p in paths[i]) for n in Destination[i]) == - orderOriginalStartTime[i] )

	#Truck constraints
	@constraint(rmp, initialTrucks[n in N_0], sum(sum(sum(delta[i,a,p] * x[i,p] for a in A_plus_i_full[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in A_plus_hd[n]) == m_0[n])
	@constraint(rmp, finalTrucks[n in N_end], sum(sum(sum(delta[i,a,p] * x[i,p] for a in A_minus_i_full[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) >= m_end[n])
	@constraint(rmp, truckFlowBalance[n in N_flow_t], sum(sum(sum(delta[i,a,p] * x[i,p] for a in A_minus_i_full[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) - sum(sum(sum(delta[i,a,p] * x[i,p] for a in A_plus_i_full[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) - sum(y[a] for a in A_plus_hd[n]) == 0)
	
	#Linking constraints
	@constraint(rmp, driverAvailability[a in A_space], sum(sum(sum(z[l,s,f] for f in fragmentscontaining[l,s,a]) for s in 1:numeffshifts) for l in 1:numlocs) == w[a] )
	for i in orders, p in paths[i], a in orderArcSet_space_full[i]
		set_normalized_coefficient(driverAvailability[a], x[i,p], -1 * delta[i,a,p])
	end
	for a in A_hasdriver_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(rmp, driverStartingLocs[l in 1:numlocs, s in 1:numeffshifts], sum(sum(z[l,s,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == length(driversets[l,s]))
	@constraint(rmp, driverFlowBalance[l in 1:numlocs, s in 1:numeffshifts, n in N_flow_ls[l,s]], sum(z[l,s,f] for f in F_minus_ls[l,s,n]) - sum(z[l,s,f] for f in F_plus_ls[l,s,n]) == 0)

	#optimize!(rmp)
	#objective_value(rmp)

	#sum(sum(sum(z[l,s,f] for f in fragmentscontaining[l,s,a]) 

	#for s in 1:numeffshifts, l in 1:numlocs, a in 1:numarcs
	#	for f in fragmentscontaining[l,s,a]
	#		println("$l, $s, $a, $f")
	#	end
	#end

	#Return named tuple of constraints needed for column generation
	rmpconstraints = (con_orderpath = orderpath,
		con_initialTrucks = initialTrucks,
		con_finalTrucks = finalTrucks,
		con_truckFlowBalance = truckFlowBalance,
		con_driverAvailability = driverAvailability,
		con_deliveryTime = deliveryTime
		)

	return rmp, x, y, z, w, rmpconstraints

end

#----------------------------------------------------------------------------------------#

function findarcvariablereducedcosts_pbcg(z, alpha, beta, epsilon, gamma, eta, psi, orderArcSet_space)

	arcredcosts = Dict()
	for i in orders, a in 1:extendednumarcs
		arcredcosts[i,a] = c[a] 
	end

	#Calculate reduced costs for each arc
	for i in orders, a in orderArcSet_space[i]
		arcredcosts[i,a] += eta[a] 
	end
	for i in orders, n in N_0, a in setdiff(A_plus_i_full[i,n], dummyarc)
		arcredcosts[i,a] -= beta[n]
	end
	for i in orders, n in N_end, a in setdiff(A_minus_i_full[i,n], dummyarc)
		arcredcosts[i,a] -= epsilon[n]
	end
	for i in orders, n in N_flow_t, a in setdiff(A_minus_i_full[i,n], dummyarc)
		arcredcosts[i,a] -= gamma[n]
	end
	for i in orders, n in N_flow_t, a in setdiff(A_plus_i_full[i,n], dummyarc)
		arcredcosts[i,a] += gamma[n]
	end
	for i in orders, n in Destination[i], a in setdiff(A_minus_i_full[i,n], dummyarc)
		#arcredcosts[i,a] += nodesLookup[arcLookup[a][2]][2] * psi[i]
		arcredcosts[i,a] += arcfinishtime[a] * psi[i]
	end

	return arcredcosts

end

#----------------------------------------------------------------------------------------#


#lambda * sum((value(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) +

#for i in orders
#	println(value(ordtime[i]))
#end
#sum(sum(sum(c[a]*delta[i,a,p]*x[i,p] for a in orderArcSet[i] ) for p in paths[i]) for i in orders) + sum(c[a]*y[a] for a in A_hasdriver) + sum(u[a]*w[a] for a in A_space) )


#for i in orders, p in intersect(paths[i], [3])
#	for a in orderArcSet[i]
#		if delta[i,a,p] > 0.01
#			println("delta[$i,$a,$p] = ", delta[i,a,p])
#		end
#	end
#end

function fragmentcolumngeneration!(paths, delta, orderArcSet_full, orderArcSet_space_full, numeffshifts)

	#Build sparse master problem
	model, x, y, z, w, rmpconstraints = sparsemasterproblem_pbcg(paths, delta, orderArcSet_full, orderArcSet_space_full, numeffshifts)

	#------------------------------------------------------#

	#Initialize column generation 
	cg_iter = 1
	rmpobjectives, rmptimes, pptimes, pptimes_par, lowerbounds = [], [], [], [], []
	listlength = convert(Int64, ceil(length(orders)/4))
	shuffle_partition(N; chunk_size=listlength) = (collect ∘ partition)(shuffle(1:N), chunk_size)

	#------------------------------------------------------#

	#Pre-processing - constraint matrix + subproblem
	#xcoeffs_alpha, xcoeffs_beta, xcoeffs_gamma, xcoeffs_theta, xcoeffs_nu, xcoeffs_mu, xcoeffs_xi, xcoeffs_psi = findxcoefficientmatrix(A_minus_i_full, A_plus_i_full) 
	numnodesSP, arclistSP, nodesLookupSP, arcLookupSP, dummyorig, dummydest = mvgsubproblem_preprocess(orderArcSet_full)

	#------------------------------------------------------#

	while cg_iter <= 100000

		#-----------SOLVE rmp-----------#

		#println("-------- ITERATION $cg_iter --------")

		status = optimize!(model)
		if termination_status(model) != MOI.OPTIMAL
			println(termination_status(model))
			return 100000000, model, x, y, z, w, paths, delta
		end
		rmpobj, rmptime = objective_value(model), solve_time(model)
		
		push!(rmpobjectives, copy(rmpobj))
		push!(rmptimes, copy(rmptime))	
		if printstatements == 1
			println("Solved rmp with objective = ", rmpobj, " in iteration $cg_iter (", sum(length(paths[i]) for i in orders), " paths)")
			#println("rmp time = ", rmptime)
		end

		#------------SNAPSHOT-----------#
		#Snapshot of current paths, before adding new paths - used for reference in warm start and column management
		paths_snapshot = deepcopy(paths)

		#Initialize list of new paths (to be added to RMP)
		newpaths = []

		#---------RED COST CALC---------#

		starttime2 = time()
	
		#Get constraint dual values
		alpha = dual.(rmpconstraints.con_orderpath)
		beta = dual.(rmpconstraints.con_initialTrucks)
		epsilon = dual.(rmpconstraints.con_finalTrucks)
		gamma = dual.(rmpconstraints.con_truckFlowBalance)
		eta = dual.(rmpconstraints.con_driverAvailability)
		psi = dual.(rmpconstraints.con_deliveryTime)

		arcredcosts = findarcvariablereducedcosts_pbcg(z, alpha, beta, epsilon, gamma, eta, psi, orderArcSet_space_full)

		#------------SUBPROBLEM------------#

		#Run shortest path for each order to find new arcs
		shortestpathnodes, shortestpatharcs = Dict(), Dict()
		dptimelist, dptimelist_ber = [], []
		addarcs, minreducedcosts = [], []
		for i in orders
			starttime = time()
			pathredcost, spnodes, sparcs = mvgsubproblem_fast(i, arcredcosts, orderArcSet_full, numnodesSP, arclistSP[i], nodesLookupSP, arcLookupSP[i], dummyorig, dummydest)
			endtime = time() - starttime

			minreducedcost = pathredcost - alpha[i]

			push!(dptimelist, endtime)
			push!(minreducedcosts, minreducedcost)

			if minreducedcost < -0.0001
				newpathindex = last(paths[i]) + 1
				pathobjcost = 0
				for a in orderArcSet_full[i]
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
		push!(pptimes, sum(dptimelist))

		#"Parallelize"
		shuffleddptimes = shuffle_partition(length(orders))
		dptimelistsums = [sum(dptimelist[j_inner] for j_inner in shuffleddptimes[j_outer]) for j_outer in 1:length(shuffleddptimes)]
		push!(pptimes_par, maximum(dptimelistsums))
		try
			#push!(lowerbounds, rmpobj + sum(minreducedcosts[k] for k in 1:length(minreducedcosts) if minreducedcosts[k] <= -0.000001))
			push!(lowerbounds, rmpobj + sum(minreducedcosts[k] for k in 1:length(minreducedcosts) if minreducedcosts[k] <= -0.000001))
		catch
			push!(lowerbounds, rmpobj)
		end
		#println("PP time (parallel)= ", maximum(dptimelistsums))

		#------COUNT ARCS AND PATHS-----#

		if saveconvergencedata_flag >= 0
			totalorderarcs = sum(sum(min(1, sum(delta[i,a,p] for p in paths_snapshot[i])) for a in orderArcSet_full[i]) for i in orders)
			totalorderpaths = sum(length(paths_snapshot[i]) for i in orders)
			maximprove = minimum(minreducedcosts) * sum(sum(value(x[i,p]) for p in paths_snapshot[i]) for i in orders)
			write_cg_conv(convergencedatafilename, cg_iter, maximprove, totalorderarcs, totalorderpaths, rmpobj)
		end
	
		#-------ADD NEW VARIABLES-------#

		for row in newpaths
			#Get order index, path index, and path cost
			i, p = row[1], row[2]
			cost = sum(c[a]*delta[i,a,p] for a in orderArcSet_full[i]) #row[3]
		
			#Create a new variable for the path
			global x[i,p] = @variable(model, upper_bound = 1, lower_bound = 0)
			set_name(x[i,p], string("x[",i,",",p,"]")) 

			#Add to objective
			set_objective_function(model, objective_function(model) + cost*x[i,p]) 

			#Add new variable to order path constraints
			set_normalized_coefficient(rmpconstraints.con_orderpath[i], x[i,p], 1.0)

			#Add to delivery time constraintsinclude("pathbasedcolumngeneration.jl")
			set_normalized_coefficient(rmpconstraints.con_deliveryTime[i], x[i,p], - sum(sum(arcfinishtime[a] * delta[i,a,p] for a in A_minus_i_full[i,n]) for n in Destination[i] if A_minus_i_full[i,n] != []) )
			
			#Add new variable to initial and final trucks constraints
			for l in 1:numlocs
				if A_plus_i_full[i,nodes[l,0]] != []
					set_normalized_coefficient(rmpconstraints.con_initialTrucks[nodes[l,0]], x[i,p], sum(delta[i,a,p] for a in A_plus_i_full[i,nodes[l,0]]))
				end
				if A_minus_i_full[i,nodes[l,horizon]] != []
					set_normalized_coefficient(rmpconstraints.con_finalTrucks[nodes[l,horizon]], x[i,p], sum(delta[i,a,p] for a in A_minus_i_full[i,nodes[l,horizon]]))
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
			for a in intersect(orderArcSet_space_full[i], shortestpatharcs[i])
				#println("Old")
				#println(rmpconstraints.con_driverAvailability[a])
				set_normalized_coefficient(rmpconstraints.con_driverAvailability[a], x[i,p], -1.0)
				#println("New")
				#println(rmpconstraints.con_driverAvailability[a])
			end

			#for a in orderArcSet_space[i]
			#	println("Old")
			#	println(rmpconstraints.con_driverAvailability[a])
			#	set_normalized_coefficient(rmpconstraints.con_driverAvailability[a], x[i,p], -1 * delta[i,a,p])
			#	println("New")
			#	println(rmpconstraints.con_driverAvailability[a])
			#end

		end	

		#----------TERMINATION----------#

		if minimum(minreducedcosts) >= -0.0001
			if printstatements == 1
				println("NO NEGATIVE REDUCED COSTS FOUND!")	
			end

			break

		end

		#------------ITERATE------------#

		cg_iter += 1

	end

	optimize!(model)
	rmpobj = objective_value(model)

	#-------------------------------#

	df = DataFrame(experiment_id = [experiment_id for j in 1:length(rmptimes)],
			instance = [ex for j in 1:length(rmptimes)],
			lambda = [lambda for j in 1:length(rmptimes)],
			horizon = [horizon for j in 1:length(rmptimes)],
			tstep = [tstep for j in 1:length(rmptimes)],
			week = [weekstart for j in 1:length(rmptimes)],
			driverfactor = [driverfactor for j in 1:length(rmptimes)],
			numdrivers = [length(drivers) for j in 1:length(rmptimes)],
			method = [solutionmethod for j in 1:length(rmptimes)],
			k = [k for j in 1:length(rmptimes)],
			iteration = [j for j in 1:length(rmptimes)],
            objective = rmpobjectives,
			milesobj = [0 for j in 1:length(rmptimes)],
			delayobj = [0 for j in 1:length(rmptimes)],
            bound = lowerbounds,
            rmptime = rmptimes,
            pptime = pptimes,
            pptime_par = pptimes_par,
            iptime = [0 for j in 1:length(rmptimes)],
			mvgtime_full = [0 for j in 1:length(rmptimes)],
           )

	CSV.write(resultsfilename, df)

	#-------------------------------#

	return rmpobj, model, x, y, z, w, paths, delta, sum(rmptimes), sum(pptimes), sum(pptimes_par)

end

#----------------------------------------------------------------------------------------#

function solvefragmentip_paths(opt_gap, paths, delta, orderArcSet, orderArcSet_space, numeffshifts)

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60*90)
	set_optimizer_attribute(ip, "OutputFlag", 1)
	set_optimizer_attribute(ip, "MIPGap", opt_gap)

	#Variables
	@variable(ip, x[i in orders, p in paths[i]], Bin)
	@variable(ip, y[A_hasdriver] >= 0, Int)
	@variable(ip, z[l = 1:numlocs, s = 1:numeffshifts, f = 1:numfragments[l,s]] >= 0, Int)
	@variable(ip, w[a in A_space] >= 0)
	@variable(ip, ordtime[orders])

	#Objective
	@objective(ip, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(sum(c[a]*delta[i,a,p]*x[i,p] for a in orderArcSet[i] ) for p in paths[i]) for i in orders) + sum(c[a]*y[a] for a in A_hasdriver) + sum(u[a]*w[a] for a in A_space) )

	#Order constraints
	@constraint(ip, orderpath[i in orders], sum(x[i,p] for p in paths[i]) == 1)

	#Order delivery constraints
	@constraint(ip, deliveryTime[i in orders], ordtime[i] - sum(sum(sum(arcfinishtime[a] * delta[i,a,p] * x[i,p] for a in A_minus_i_full[i,n]) for p in paths[i]) for n in Destination[i]) == - orderOriginalStartTime[i] )

	#Truck constraints
	@constraint(ip, initialTrucks[n in N_0], sum(sum(sum(delta[i,a,p] * x[i,p] for a in A_plus_i_full[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in A_plus_hd[n]) == m_0[n])
	@constraint(ip, finalTrucks[n in N_end], sum(sum(sum(delta[i,a,p] * x[i,p] for a in A_minus_i_full[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) >= m_end[n])
	@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(sum(delta[i,a,p] * x[i,p] for a in A_minus_i_full[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) - sum(sum(sum(delta[i,a,p] * x[i,p] for a in A_plus_i_full[i,n]) for p in setdiff(paths[i], dummypath)) for i in orders) - sum(y[a] for a in A_plus_hd[n]) == 0)
	
	#Linking constraints
	@constraint(ip, driverAvailability[a in A_space], sum(sum(sum(z[l,s,f] for f in fragmentscontaining[l,s,a]) for s in 1:numeffshifts) for l in 1:numlocs) == w[a] )
	for i in orders, p in paths[i], a in orderArcSet_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,p], -1 * delta[i,a,p])
	end
	for a in A_hasdriver_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(ip, driverStartingLocs[l in 1:numlocs, s in 1:numeffshifts], sum(sum(z[l,s,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == length(driversets[l,s]))
	@constraint(ip, driverFlowBalance[l in 1:numlocs, s in 1:numeffshifts, n in N_flow_ls[l,s]], sum(z[l,s,f] for f in F_minus_ls[l,s,n]) - sum(z[l,s,f] for f in F_plus_ls[l,s,n]) == 0)

	optimize!(ip)

	ip_obj = objective_value(ip)
	println("Frag path-based IP objective = ", ip_obj)
	println("Time = ", solve_time(ip))

	return ip_obj, value.(z), solve_time(ip)

end

