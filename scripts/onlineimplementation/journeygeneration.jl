
function getdualvalues_journey(rmpconstraints)

	alph = Array(dual.(rmpconstraints.con_driverAvailability))
	beta = Array(dual.(rmpconstraints.con_driverStartingLocs))
	gamma = dual.(rmpconstraints.con_driverFlowBalance)

    return alph, beta, gamma

end

#=
reduced_cost.(z)
(hl,ss,sn,lth), j = (17, 2, 17, 0), 9 

reduced_cost(z[(hl,ss,sn,lth), j]) = 988.7936717205016
myarcs = currfragments.fragmentarcs[hl,ss,sn,lth, j]
currfragments.fragmentarcs[hl,ss,sn,lth, j] = [355, 6748, 6749, 6750, 977, 7056, 7057, 7058, 5792]
myarcs = [355, 6748, 6749, 6750, 977, 7056, 7057, 7058, 5792]
n_start = arcLookup[myarcs[1]][1]
n_end = arcLookup[last(myarcs)][2]
g = currfragments.drivergroupnum[hl,ss,sn,lth]

beta[g] + sum(arcredcosts[1,a] for a in myarcs) - gamma[(hl,ss,sn,lth), n_end] + gamma[(hl,ss,sn,lth), n_start]
for a in myarcs
	println("$a --> ", arcredcosts[a])
	arcDesc(a)
end
=#

#----------------------------------------------------------------------------------------#

function findjourneyarcreducedcosts(rmpconstraints, currfragments, spacevec)

	alph, beta, gamma = getdualvalues_journey(rmpconstraints)
	
	mapdriversettoindex = Dict()
	for i in 1:length(currfragments.driversets)
		mapdriversettoindex[currfragments.driversets[i]] = i
	end

	arcredcosts = zeros(length(currfragments.driversets), extendednumarcs)
	for (hl,ss,sn,lth) in currfragments.driversets
		arcredcosts[mapdriversettoindex[hl,ss,sn,lth],:] -= spacevec * alph
	end

	return arcredcosts, beta, gamma, mapdriversettoindex

end

#----------------------------------------------------------------------------------------#

function findinitialjourneyarcs(currstate, hl, ss, sn, lth, d_ex, ghostdriverarcs, offtimestarttimes)

	#Find deadline and round down to beginning of an off period
	if (sn == nodes[hl,0]) & (0 in currstate.T_off[ss])
		deadlinetime = minimum(currstate.T_on_0[d_ex])
	else
		deadlinetime = lth + (shiftlength + 24*maxnightsaway)
	end
	#offtimestarttimes = currstate.T_off_0[d_ex] #union(currstate.T_off_0[d_ex], last(currstate.T_off_0[d_ex])+24:24:deadlinetime)
	deadlinetime = maximum([t for t in offtimestarttimes if t <= deadlinetime]) + (24-shiftlength)
	if deadlinetime > horizon
		deadlinenode = extendednodes[hl, dummyendtime]
	else
		deadlinenode = extendednodes[hl, deadlinetime]
	end

	#Reduce arcs for the driver group based on the deadline
	availablearcs_ghost = shortestpathonTSN(sn, hl, deadlinetime, ghostdriverarcs.A[hl,ss], ghosttsn)
	availablearcs_base = translatearcsbetweenTSNs_withextendedarcs(availablearcs_ghost, ghosttsn, basetsn)
	journeyarcs = formatarcset(availablearcs_base, basetsn)

	return journeyarcs

end

#----------------------------------------------------------------------------------------#

function findfuturejourneyarcs(t, currstate, hl, ss, sn, lth, ghostdriverarcs, offtimestarttimes)

	fsn = nodes[hl, t] #Future start node
	future_deadlinetime = t + shiftlength + 24*maxnightsaway
	future_deadlinetime = maximum([t for t in offtimestarttimes if t <= future_deadlinetime]) + (24-shiftlength)
	if future_deadlinetime > horizon
		future_deadlinenode = extendednodes[hl, dummyendtime]
	else
		future_deadlinenode = extendednodes[hl, future_deadlinetime]
	end

	future_availablearcs_ghost = shortestpathonTSN(fsn, hl, future_deadlinetime, ghostdriverarcs.A[hl,ss], ghosttsn)
	future_availablearcs_base = translatearcsbetweenTSNs_withextendedarcs(future_availablearcs_ghost, ghosttsn, basetsn)
	
	future_journeyarcs = formatarcset(future_availablearcs_base, basetsn)

	return future_journeyarcs

end

#----------------------------------------------------------------------------------------#

function findalljourneyarcs(currstate, currfragments, hl, ss, sn, lth, ghostdriverarcs)

	journeyarcsfor = Dict()

	startloc, starttime = nodesLookup[sn]
	g = currfragments.drivergroupnum[hl,ss,sn,lth]
	d_ex = currfragments.driversingroup[hl,ss,sn,lth][1] #Example driver for all the annoying sets that are indexed by driver instead of group
	upcomingfragmentstarttime = [t for t in currstate.T_on_0[d_ex] if (t > starttime) & (t < horizon)]
	offtimestarttimes = currstate.T_off_0[d_ex]
	
	journeyarcs_init = findinitialjourneyarcs(currstate, hl, ss, sn, lth, d_ex, ghostdriverarcs, offtimestarttimes)
	journeyarcsfor[basetsn.nodesLookup[sn][2]] = journeyarcs_init

	for t in upcomingfragmentstarttime
		journeyarcs = findfuturejourneyarcs(t, currstate, hl, ss, sn, lth, ghostdriverarcs, offtimestarttimes)
		journeyarcsfor[t] = journeyarcs
	end

	return journeyarcsfor

end

#----------------------------------------------------------------------------------------#

function journeyshortestpathproblem(arccost, spsets)

	startnode = spsets.startnode
	
	#Initialize shortest path algorithm 
	currdistance = repeat([999999999.0],outer=[spsets.numnodes])
	currdistance[startnode] = 0
	prevnode, prevarc = zeros(spsets.numnodes), zeros(spsets.numnodes)
	
	#Loop over time-space arcs in order of start time
	for a in spsets.availablearcs
		#println(a)
		n_end, n_start = spsets.arclookup[a][2], spsets.arclookup[a][1]
		if currdistance[n_end] > currdistance[n_start] + arccost[a] + 0.000001
			currdistance[n_end] = currdistance[n_start] + arccost[a]
			prevnode[n_end] = n_start
			prevarc[n_end] = a
		end
	end

	#setdiff(spsets.availablearcs, 1:extendednumarcs)
	#for a in setdiff(spsets.availablearcs, 1:extendednumarcs)
	#	println("$a --> ", spsets.arclookup[a])
	#end
	

	#myarcs = [a for a in intersect(spsets.availablearcs, 1:numarcs) if currdistance[arcLookup[a][1]] <= 100000]
	#timespacenetwork("outputs/viz/aaa_all.png", [myarcs], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
		

	#Format the shortest path output
	shortestpathnodes_rev = [spsets.dummydest]
	shortestpatharcs_rev = []
	node = spsets.dummydest
	while node != startnode
		push!(shortestpatharcs_rev, Int(prevarc[node]))
		node = Int(prevnode[node])
		push!(shortestpathnodes_rev, node)
	end
	shortestpathnodes = reverse(shortestpathnodes_rev[2:length(shortestpathnodes_rev)]) 
	shortestpatharcs = reverse(shortestpatharcs_rev[2:length(shortestpatharcs_rev)]) 

	return currdistance[spsets.dummydest], shortestpathnodes, shortestpatharcs

end

#----------------------------------------------------------------------------------------#

#hl,ss,sn,lth, arcredcosts, beta, gamma, journeyarcsfor, currfragments, subproblemsets = hl,ss,sn,lth, arcredcosts, beta, gamma, alljourneyarcs[hl,ss,sn,lth], currfragments, subproblemsets

function journeysubproblem(hl,ss,sn,lth, arcredcosts, beta, gamma, journeyarcsfor, currfragments, subproblemsets)	

	newjourneys, journeyreducedcosts = [], []
	minrc = 1e10

	journeystarttimes = keys(journeyarcsfor)
	g = currfragments.drivergroupnum[hl,ss,sn,lth]

	for t in journeystarttimes
		#Get start node
		#if t == nodesLookup[sn][2]
		#	startcost = beta[g]
		#else
		#	startcost = 0
		#end

		#Finalize arc costs
		arccost = Dict()
		for a in subproblemsets[hl,ss,sn,lth,t].availablearcs
			if a <= extendednumarcs
				arccost[a] = arcredcosts[g,a]
			else
				arccost[a] = 0
			end
		end
		for a in subproblemsets[hl,ss,sn,lth,t].A_minus_dummydest
			n = subproblemsets[hl,ss,sn,lth,t].arclookup[a][1]
			if n in currfragments.N_flow_g[hl,ss,sn,lth]
				arccost[a] -= gamma[(hl,ss,sn,lth), n]
			end
		end

		#Find the min reduced cost journey 
		journeyreducedcost, journeynodelist, journeyarclist = journeyshortestpathproblem(arccost, subproblemsets[hl,ss,sn,lth,t])
		if t == nodesLookup[sn][2]
			journeyreducedcost = journeyreducedcost + beta[g] + 0
		else
			journeyreducedcost = journeyreducedcost + 0 + gamma[(hl,ss,sn,lth), subproblemsets[hl,ss,sn,lth,t].startnode]
		end

		#Add promising journeys
		if journeyreducedcost < -1e-4
			push!(newjourneys, journeyarclist)
			push!(journeyreducedcosts, journeyreducedcost)
		end

		minrc = min(minrc, journeyreducedcost)
				
	end

	return newjourneys, journeyreducedcosts, minrc

end

#----------------------------------------------------------------------------------------#

function preprocesssubproblemsets_journeygen(alljourneyarcs, currfragments)

	subproblemsets = Dict()

	for (hl,ss,sn,lth) in currfragments.driversets
		journeyarcsfor = alljourneyarcs[hl,ss,sn,lth]
		journeystarttimes = keys(journeyarcsfor)
		d_ex = currfragments.driversingroup[hl,ss,sn,lth][1] 

		for t in journeystarttimes
			#Get start node
			if t == nodesLookup[sn][2]
				journeystartnode = sn
			else
				journeystartnode = nodes[hl,t]
			end

			#Finalize arc costs
			sp_availablearcs = journeyarcsfor[t].A
			journeyendingnodelist = union(N_end, [nodes[hl, t] for t in [t2 for t2 in currstate.T_on_0[d_ex] if t2 <= horizon] if nodes[hl,t] > sn], numnodes+1:extendednumnodes)
			sp_dummydest = extendednumnodes + 1
			sp_numnodes = extendednumnodes + 1
			sp_arclookup = copy(arcLookup)
			A_minus_dummydest = []
			
			arcindex = dummyarc+1
			for n in journeyendingnodelist
				sp_arclookup[arcindex] = (n,sp_dummydest)
				push!(sp_availablearcs, arcindex)
				push!(A_minus_dummydest, arcindex)
				arcindex += 1
			end

			subproblemsets[hl,ss,sn,lth,t] = (numnodes=sp_numnodes, arclookup=sp_arclookup, dummydest=sp_dummydest, startnode=journeystartnode, availablearcs=sp_availablearcs, A_minus_dummydest=A_minus_dummydest)
			
		end
	end

	return subproblemsets

end

#----------------------------------------------------------------------------------------#

#opt_gap, orderarcs, ghostdriverarcs = opt_gap, currarcs.orderarcs, currarcs.ghostdriverarcs

function journeygeneration(opt_gap, orderarcs, ghostdriverarcs)
	
    #include("scripts/onlineimplementation/journeygeneration.jl")
	#include("scripts/onlineimplementation/initializecurrentstatearcs.jl")
	#currarcs, currfragments, primaryarcs, extendedtimearcs, numeffshifts, journeytime = initializecurrentstatearcs(currstate, 0 )

	currjourneys = Dict()
    for (hl,ss,sn,lth) in currfragments.driversets
        currjourneys[hl,ss,sn,lth] = [j for j in 1:currfragments.numfragments[hl,ss,sn,lth]]
    end

    #-------------------------------------------------------------------------------#

	rmp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(rmp, "TimeLimit", 60*60*40)
	set_optimizer_attribute(rmp, "OutputFlag", 0)
	set_optimizer_attribute(rmp, "MIPGap", opt_gap)

	A_space_all = primaryarcs.A_space
	for i in currstate.orders
		goodones = [a for a in orderarcs.A_space[i] if nodesLookup[arcLookup[a][1]][2] < horizon]
		A_space_all = union(A_space_all, goodones)
	end

	#Variables
	@variable(rmp, 0 <= x[i in currstate.orders, a in orderarcs.A[i]] )
	@variable(rmp, y[currarcs.hasdriverarcs.A] >= 0)
	@variable(rmp, z[(i1,i2,i3,i4) = currfragments.driversets, f = currjourneys[i1,i2,i3,i4]] >= 0)
	@variable(rmp, w[a in A_space_all] >= 0)
	@variable(rmp, ordtime[currstate.orders])

	#Objective
	@objective(rmp, Min, lambda * sum((ordtime[i] - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i] for i in currstate.orders) 
	+ sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in currstate.orders) + sum(c[a]*(y[a] ) for a in currarcs.hasdriverarcs.A) + sum(u[a]*(w[a] ) for a in primaryarcs.A_space) )

	#Order constraints
	@constraint(rmp, orderFlowBalance[i = currstate.orders, n in setdiff([n2 for n2 in 1:numnodes], union(currstate.Origin[i], currstate.Destination[i]))], sum(x[i,a] for a in orderarcs.A_minus[i,n]) - sum(x[i,a] for a in orderarcs.A_plus[i,n]) == 0)
	@constraint(rmp, arriveDestin[i = currstate.orders], sum(sum(x[i,a] for a in orderarcs.A_minus[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in currstate.Destination[i]) == 1)
	@constraint(rmp, departOrigin[i = currstate.orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), orderarcs.A_plus[i,n])) for n in currstate.Origin[i]) == 1)
	for i in setdiff(currstate.orders, currstate.ordersinprogress)
		extendedorderarc = extendedarcs[last(currstate.Origin[i]), last(currstate.Destination[i])]
		if extendedorderarc in orderarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end

	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(currstate.orders, currstate.ordersinprogress), n in currstate.Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
		if a in orderarcs.A[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(rmp, deliveryTime[i in currstate.orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in orderarcs.A_minus[i,n]) for n in currstate.Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(rmp, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == currstate.m_0[n])
	@constraint(rmp, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) >= currstate.m_end[n])
	#@constraint(rmp, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n],dummyarc)) for i in currstate.orders) - sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == 0)
	@constraint(rmp, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in currstate.orders) + sum(y[a] for a in currarcs.hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n],dummyarc)) for i in currstate.orders) - sum(y[a] for a in currarcs.hasdriverarcs.A_plus[n]) == 0)

	#Linking constraints
	@constraint(rmp, driverAvailability[a in A_space_all], sum(sum(z[(i1,i2,i3,i4),f] for f in currfragments.fragmentscontaining[i1,i2,i3,i4,a]) for (i1,i2,i3,i4) in currfragments.driversets) == w[a]  )
	for i in currstate.orders, a in [a for a in orderarcs.A_space[i] if nodesLookup[arcLookup[a][1]][2] < horizon]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in currarcs.hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(rmp, driverStartingLocs[(i1,i2,i3,i4) in currfragments.driversets], sum(sum(z[(i1,i2,i3,i4),f] for f in currfragments.F_plus_g[i1,i2,i3,i4,n]) for n in [i3]) == length(currfragments.driversingroup[i1,i2,i3,i4]))
	@constraint(rmp, driverFlowBalance[(i1,i2,i3,i4) in currfragments.driversets, n in currfragments.N_flow_g[i1,i2,i3,i4]], sum(z[(i1,i2,i3,i4),f] for f in currfragments.F_minus_g[i1,i2,i3,i4,n]) - sum(z[(i1,i2,i3,i4),f] for f in currfragments.F_plus_g[i1,i2,i3,i4,n]) == 0)

	#Return named tuple of constraints needed for column generation
	rmpconstraints = (
		con_driverAvailability = driverAvailability,
		con_driverStartingLocs = driverStartingLocs,
		con_driverFlowBalance = driverFlowBalance
		)

	#-------------------------------------------------------------------------------#

	spacevec = zeros(extendednumarcs,length(A_space_all))
	for aindex in 1:length(A_space_all)
		spacevec[A_space_all[aindex],aindex] = 1
	end

	rmpobjectives, rmptimes = [], []

	include("scripts/onlineimplementation/journeygeneration.jl")

	alljourneyarcs = Dict()
	for (hl,ss,sn,lth) in currfragments.driversets
		journeyarcsfor = findalljourneyarcs(currstate, currfragments, hl, ss, sn, lth, ghostdriverarcs)
		alljourneyarcs[hl,ss,sn,lth] = journeyarcsfor
	end

	subproblemsets = preprocesssubproblemsets_journeygen(alljourneyarcs, currfragments)

	#-------------------------------------------------------------------------------#

	cg_iter = 1

	while cg_iter <= 100000
				
		#-------------SOLVE RMP-------------#
		println("-------- ITERATION $cg_iter --------")
		@time status = optimize!(rmp)
		if termination_status(rmp) != MOI.OPTIMAL
			println(termination_status(rmp))
			return 100000000, rmp, x, y, z, w, paths, delta
		end
		rmpobj, rmptime = objective_value(rmp), solve_time(rmp)
		push!(rmpobjectives, copy(rmpobj))
		push!(rmptimes, copy(rmptime))	
		println("Solved RMP with objective = ", rmpobj, " in iteration $cg_iter (", sum(length(currjourneys[hl,ss,sn,lth]) for (hl,ss,sn,lth) in currfragments.driversets), " journeys)")

		#------------SUBPROBLEM-------------#

		#Calculate reduced costs
		arcredcosts, beta, gamma, mapdriversettoindex = findjourneyarcreducedcosts(rmpconstraints, currfragments, spacevec)

		#Run shortest path for each order to find new arcs
		shortestpathnodes, shortestpatharcs = Dict(), Dict()
		dptimelist, minreducedcosts = [], []
		for (hl,ss,sn,lth) in currfragments.driversets

			newjourneys, journeyreducedcosts, minreducedcost = journeysubproblem(hl,ss,sn,lth, arcredcosts, beta, gamma, alljourneyarcs[hl,ss,sn,lth], currfragments, subproblemsets)

			#timespacenetwork("outputs/viz/aaa_all.png", [alljourneyarcs[hl,ss,sn,lth][72].A], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
			#timespacenetwork("outputs/viz/aaa_all.png", [newjourneys[1]], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)	

			push!(minreducedcosts, minreducedcost)

			for journeryarclist in newjourneys

				#Add new journey
				j = last(currjourneys[hl,ss,sn,lth]) + 1
				push!(currjourneys[hl,ss,sn,lth], j)

				#Add to sets
				n_1,n_2 = arcLookup[journeryarclist[1]][1], arcLookup[last(journeryarclist)][2]
				push!(currfragments.F_minus_g[hl,ss,sn,lth,n_2], j)
				push!(currfragments.F_plus_g[hl,ss,sn,lth,n_1], j)
				currfragments.fragmentarcs[hl,ss,sn,lth,j] = journeryarclist
				for a in journeryarclist
					push!(currfragments.fragmentscontaining[hl,ss,sn,lth,a], j)
				end

				#-------ADD NEW VARIABLES-------#

				#Create a new variable for the path
				global z[(hl,ss,sn,lth),j] = @variable(rmp, lower_bound = 0)
				set_name(z[(hl,ss,sn,lth),j], string("z[(",hl,",",ss,",",sn,",",lth,"),",j,"]")) 
			
				#Linking constraints
				for a in intersect(journeryarclist, A_space_all)
					set_normalized_coefficient(rmpconstraints.con_driverAvailability[a], z[(hl,ss,sn,lth),j], 1.0)
				end
			
				#Driver constraints
				#println(rmpconstraints.con_driverStartingLocs)
				if sn == n_1
					set_normalized_coefficient(rmpconstraints.con_driverStartingLocs[(hl,ss,sn,lth)], z[(hl,ss,sn,lth),j], 1.0)
				end
				if n_1 in currfragments.N_flow_g[hl,ss,sn,lth]
					set_normalized_coefficient(rmpconstraints.con_driverFlowBalance[(hl,ss,sn,lth),n_1], z[(hl,ss,sn,lth),j], -1.0)
				end
				if n_2 in currfragments.N_flow_g[hl,ss,sn,lth]
					set_normalized_coefficient(rmpconstraints.con_driverFlowBalance[(hl,ss,sn,lth),n_2], z[(hl,ss,sn,lth),j], 1.0)
				end
			
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

	return objective_value(rmp), value.(x), value.(z), value.(w), value.(y)
		
end




#=

include("scripts/onlineimplementation/journeygeneration.jl")
currarcs, currfragments, primaryarcs, extendedtimearcs, numeffshifts, journeytime = initializecurrentstatearcs(currstate, 0);
lp_obj, x_lp, z_lp, w_lp, y_lp = journeygeneration(opt_gap, currarcs.orderarcs, currarcs.ghostdriverarcs)

=#