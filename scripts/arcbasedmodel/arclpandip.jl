
function solvelp(orderArcSet, orderArcSet_space, A_plus_i, A_minus_i)

	lp = Model(Gurobi.Optimizer)
	set_optimizer_attribute(lp, "TimeLimit", 60*60*47)
	set_optimizer_attribute(lp, "OutputFlag", 0)

	#Variables
	@variable(lp, 0 <= x[i in orders, a in orderArcSet[i]] <= 1)
	@variable(lp, y[A_hasdriver] >= 0)
	@variable(lp, 0 <= z[d = drivers, homeArcSet[d]] <= 1)
	@variable(lp, w[a in A_space] >= 0)
	@variable(lp, ordtime[orders])

	#Objective
	@objective(lp, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderArcSet[i]) for i in orders) + sum(c[a]*y[a] for a in A_hasdriver) + sum(u[a]*w[a] for a in A_space) )

	#Order constraints
	@constraint(lp, orderFlowBalance[i = orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))], sum(x[i,a] for a in A_minus_i[i,n]) - sum(x[i,a] for a in A_plus_i[i,n]) == 0)
	@constraint(lp, arriveDestin[i = orders], sum(sum(x[i,a] for a in A_minus_i[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in Destination[i]) == 1)
	@constraint(lp, departOrigin[i = orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), A_plus_i[i,n])) for n in Origin[i]) == 1)
	for i in setdiff(orders, ordersinprogress)
		extendedorderarc = extendedarcs[last(Origin[i]), last(Destination[i])]
		if extendedorderarc in orderArcSet[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end
	
	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(orders, ordersinprogress), n in Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
		if a in orderArcSet[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(lp, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in A_minus_i[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(lp, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_plus_hd[n]) == m_0[n])
	@constraint(lp, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) >= m_end[n])
	@constraint(lp, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n])- sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n],dummyarc)) for i in orders) - sum(y[a] for a in A_plus_hd[n]) == 0)

	#Driver constraints
	@constraint(lp, driverAvailability[a in A_space], sum(z[d,a] for d in intersect(drivers, availableDrivers[a])) == w[a])
	for i in orders, a in orderArcSet_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in A_hasdriver_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end
	@constraint(lp, driverFlowBalance[d in drivers, n in N_flow_d[d]], sum(z[d,a] for a in A_minus_d[d,n]) - sum(z[d,a] for a in A_plus_d[d,n]) == 0)
	@constraint(lp, driverStartingLocs[d in drivers], sum(sum(z[d,a] for a in A_plus_d[d,n]) for n in driverStartNodes[d]) == 1)
	#@constraint(ip, returnHome[d in drivers, t in setdiff(T_off_constr[d], horizon-24:tstep:horizon)], sum(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]] for t2 in intersect([t3 for t3 in t:tstep:t+24],T_off_0[d])) >= 1)  
	#@constraint(ip, returnHomeEnd[d in drivers, t in intersect(T_off_constr[d], horizon-24:tstep:horizon)], sum(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]] for t2 in intersect([t3 for t3 in t:tstep:t+24],T_off_0[d])) + sum(z[d,a] for a in A_minus_d[d, nodes[driverHomeLocs[d], horizon]]) >= 1)  
	#@constraint(ip, returnHome[d in drivers, t in setdiff(T_off_constr[d], horizon-24:tstep:horizon)], sum(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]] for t2 in setdiff(T_off_0[d], horizon)) >= 1)  
	#@constraint(ip, returnHomeEnd[d in drivers, t in intersect(T_off_constr[d], horizon-24:tstep:horizon)], sum(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]] for t2 in setdiff(T_off_0[d], horizon)) + sum(z[d,a] for a in A_minus_d[d, nodes[driverHomeLocs[d], horizon]]) >= 1)  
	#@constraint(ip, returnHome_online[d in drivers], z[d,arcs[nodes[(driverHomeLocs[d],setdiff(T_off_0[d],0)[1])], nodes[(driverHomeLocs[d],setdiff(T_off_0[d],0)[1]+tstep)]]] >= awaylastnight[d])
	@constraint(lp, returnHome[d in drivers, arcset in R_off[d]], sum(z[d,a] for a in arcset) >= 1)  

	optimize!(lp)

	lp_obj = getobjectivevalue(lp)
	println("LP objective = ", lp_obj)
	println("Time = ", solve_time(lp))

	return lp_obj, getvalue.(z), getvalue.(x), getvalue.(y), getvalue.(w), solve_time(lp)

end

#----------------------------------------------------------------------------------------#

function solveip(orderArcSet, orderArcSet_space, A_plus_i, A_minus_i)

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60*47)
	set_optimizer_attribute(ip, "OutputFlag", 1)

	#Variables
	@variable(ip, x[i in orders, a in orderArcSet[i]], Bin)
	@variable(ip, y[A_hasdriver] >= 0, Int)
	@variable(ip, z[d = drivers, homeArcSet[d]], Bin)
	@variable(ip, w[a in A_space] >= 0)
	@variable(ip, ordtime[orders])

	#Objective
	@objective(ip, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderArcSet[i]) for i in orders) + sum(c[a]*y[a] for a in A_hasdriver) + sum(u[a]*w[a] for a in A_space) )

	#Order constraints
	@constraint(ip, orderFlowBalance[i = orders, n in setdiff([n2 for n2 in 1:numnodes], union(Origin[i], Destination[i]))], sum(x[i,a] for a in A_minus_i[i,n]) - sum(x[i,a] for a in A_plus_i[i,n]) == 0)
	@constraint(ip, arriveDestin[i = orders], sum(sum(x[i,a] for a in A_minus_i[i,n] if ((a == dummyarc) || (nodesLookup[arcLookup[a][1]][1] != nodesLookup[arcLookup[a][2]][1]))) for n in Destination[i]) == 1)
	@constraint(ip, departOrigin[i = orders], sum(sum(x[i,a] for a in intersect(union(A_space, dummyarc), A_plus_i[i,n])) for n in Origin[i]) == 1)
	for i in setdiff(orders, ordersinprogress)
		extendedorderarc = extendedarcs[last(Origin[i]), last(Destination[i])]
		if extendedorderarc in orderArcSet[i]
			set_normalized_coefficient(departOrigin[i], x[i,extendedorderarc], 1)
		end
	end
	
	#Add in "stay where you are" arc for each in transit order - MAINLY NEEDED FOR ONLINE IMPLEMENTATION
	for i in intersect(orders, ordersinprogress), n in Origin[i], a in setdiff(A_plus[n], union(A_space, dummyarc))
		if a in orderArcSet[i]
			set_normalized_coefficient(departOrigin[i], x[i,a], 1)
		end
	end

	#Order delivery constraints
	@constraint(ip, deliveryTime[i in orders], ordtime[i] - sum(sum(arcfinishtime[a] * x[i,a] for a in A_minus_i[i,n]) for n in Destination[i]) == - orderOriginalStartTime[i])

	#Truck constraints
	@constraint(ip, initialTrucks[n in N_0], sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_plus_hd[n]) == m_0[n])
	@constraint(ip, finalTrucks[n in N_end], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n]) >= m_end[n])
	@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(y[a] for a in A_minus_hd[n])- sum(sum(x[i,a] for a in setdiff(A_plus_i[i,n],dummyarc)) for i in orders) - sum(y[a] for a in A_plus_hd[n]) == 0)

	#Driver constraints
	@constraint(ip, driverAvailability[a in A_space], sum(z[d,a] for d in intersect(drivers, availableDrivers[a])) == w[a])
	for i in orders, a in orderArcSet_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in A_hasdriver_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end
	@constraint(ip, driverFlowBalance[d in drivers, n in N_flow_d[d]], sum(z[d,a] for a in A_minus_d[d,n]) - sum(z[d,a] for a in A_plus_d[d,n]) == 0)
	@constraint(ip, driverStartingLocs[d in drivers], sum(sum(z[d,a] for a in A_plus_d[d,n]) for n in driverStartNodes[d]) == 1)
	#@constraint(ip, returnHome[d in drivers, t in setdiff(T_off_constr[d], horizon-24:tstep:horizon)], sum(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]] for t2 in intersect([t3 for t3 in t:tstep:t+24],T_off_0[d])) >= 1)  
	#@constraint(ip, returnHomeEnd[d in drivers, t in intersect(T_off_constr[d], horizon-24:tstep:horizon)], sum(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]] for t2 in intersect([t3 for t3 in t:tstep:t+24],T_off_0[d])) + sum(z[d,a] for a in A_minus_d[d, nodes[driverHomeLocs[d], horizon]]) >= 1)  
	#@constraint(ip, returnHome[d in drivers, t in setdiff(T_off_constr[d], horizon-24:tstep:horizon)], sum(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]] for t2 in setdiff(T_off_0[d], horizon)) >= 1)  
	#@constraint(ip, returnHomeEnd[d in drivers, t in intersect(T_off_constr[d], horizon-24:tstep:horizon)], sum(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]] for t2 in setdiff(T_off_0[d], horizon)) + sum(z[d,a] for a in A_minus_d[d, nodes[driverHomeLocs[d], horizon]]) >= 1)  
	#@constraint(ip, returnHome_online[d in drivers], z[d,arcs[nodes[(driverHomeLocs[d],setdiff(T_off_0[d],0)[1])], nodes[(driverHomeLocs[d],setdiff(T_off_0[d],0)[1]+tstep)]]] >= awaylastnight[d])
	@constraint(ip, returnHome[d in drivers, arcset in R_off[d]], sum(z[d,a] for a in arcset) >= 1)  

	#if drivervalidinequalities_flag == 1
	#	for d in drivers
	#		for t in setdiff(0:tstep:horizon, T_off[drivershift[d]])
	#			possiblehomearrivals = [t2+tstep for t2 in setdiff(0:tstep:horizon, T_off[drivershift[d]]) if (t2 >= t) & (t2+tstep <= horizon) & (t2-t <= 36)]
	#			if (possiblehomearrivals != []) && (length([t3 for t3 in setdiff(T_on_0[d],horizon) if t3 > t + 0.001]) >= 1) #&& (sum(getvalue(z[d,a]) for a in intersect(A_space, A_plus_d[d, nodes[driverHomeLocs[d], t]])) > sum(getvalue(sum(z[d,a]) for a in intersect(A_space, A_minus_d[d, nodes[driverHomeLocs[d], t2]])) for t2 in possiblehomearrivals) )
	#				@constraint(ip, sum(z[d,a] for a in intersect(A_space, A_plus_d[d, nodes[driverHomeLocs[d], t]])) <= sum(sum(z[d,a] for a in intersect(A_space, A_minus_d[d, nodes[driverHomeLocs[d], t2]])) for t2 in possiblehomearrivals) )
	#			end
	#		end
	#	end
	#end  

	#=

	for (hl,ss,sn,aln) in driversets[1:15]
		for (a,num) in setarcs[hl,ss,sn,aln]
			@constraint(ip, sum(z[d,a] for d in driversingroup[(hl,ss,sn,aln)]) >= num)
		end
	end
	optimize!(ip)
	ip_obj = getobjectivevalue(ip)



	(hl,ss,sn,aln) = driversets[16]

	setarcs[hl,ss,sn,aln]
	for (a,num) in setarcs[hl,ss,sn,aln]
		arcDesc(a)
	end

	for a in myarcs
		arcDesc(a)
	end

	myarcs = sort([a for (a,num) in setarcs[hl,ss,sn,aln]], by=x->nodesLookup[arcLookup[x][1]][2])
	=#

	optimize!(ip)

	ip_obj = getobjectivevalue(ip)
	println("IP objective = ", ip_obj)
	println("Time = ", solve_time(ip))

	return ip_obj, getvalue.(z), getvalue.(x), getvalue.(y), getvalue.(w), solve_time(ip), objective_bound(ip)

end