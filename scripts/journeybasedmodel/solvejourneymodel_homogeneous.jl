
function getnonzeroarcs(x, orderarcs)

    orderArcSet_red, orderArcSet_space_red, A_plus_i_red, A_minus_i_red = Dict(), Dict(), Dict(), Dict()

    for i in orders
        orderArcSet_red[i] = [dummyarc]
        orderArcSet_space_red[i] = []
    end
    for i in orders, n in 1:extendednumnodes
        if n == Origin[i][1]
            A_plus_i_red[i,n] = [dummyarc]
            A_minus_i_red[i,n] = []
        elseif n == last(Destination[i])
            A_plus_i_red[i,n] = []
            A_minus_i_red[i,n] = [dummyarc]
        else 
            A_plus_i_red[i,n] = []
            A_minus_i_red[i,n] = []
        end
    end

    for i in orders
        destinationlocation = nodesLookup[Destination[i][1]][1]
        #push!(Destination[i], extendednodes[destinationlocation, dummyendtime])
        A_minus_i_red[i, extendednodes[destinationlocation, dummyendtime]] = []
        for n2 in N_end
            arc_ext = extendedarcs[n2, extendednodes[destinationlocation, dummyendtime]]
            push!(orderArcSet_red[i], arc_ext)
            push!(A_plus_i_red[i, n2], arc_ext)
            push!(A_minus_i_red[i, extendednodes[destinationlocation, dummyendtime]], arc_ext)
        end
    end

    for i in orders, a in orderarcs.A[i]
        if x[i,a] > 1e-4
            orderArcSet_red[i] = union(orderArcSet_red[i], a)
            if a in primaryarcs.A_space
                orderArcSet_space_red[i] = union(orderArcSet_space_red[i], a)
            end
        end
    end
    
    #Create A_plus and A_minus lists
    for i in orders, n in 1:numnodes, a in orderarcs.A_plus[i,n]
        if (a in orderArcSet_red[i]) #& !(a in A_plus_i_red[i,n])
            push!(A_plus_i_red[i,n], a)
        end
    end
    for i in orders, n in 1:numnodes, a in orderarcs.A_minus[i,n]
        if (a in orderArcSet_red[i]) #& !(a in A_minus_i_red[i,n])
            push!(A_minus_i_red[i,n], a)
        end
    end

	for i in orders
        orderArcSet_red[i] = unique(orderArcSet_red[i])
        orderArcSet_space_red[i] = unique(orderArcSet_space_red[i])
    end
    for i in orders, n in 1:extendednumnodes
    	A_plus_i_red[i,n] = unique(A_plus_i_red[i,n])
        A_minus_i_red[i,n] = unique(A_minus_i_red[i,n])
    end

    return orderArcSet_red, orderArcSet_space_red, A_plus_i_red, A_minus_i_red

end

#---------------------------------------------------------------------------------------#

function solvejourneymodel(lprelax_flag, opt_gap, orderarcs, numeffshifts, cuts)

	ip = Model(Gurobi.Optimizer)
	set_optimizer_attribute(ip, "TimeLimit", 60*60*40)
	set_optimizer_attribute(ip, "OutputFlag", 1)
	set_optimizer_attribute(ip, "MIPGap", opt_gap)

	#Variables
    if lprelax_flag == 0
        @variable(ip, x[i in orders, a in orderarcs.A[i]], Bin)
        @variable(ip, y[hasdriverarcs.A] >= 0, Int)
        @variable(ip, z[l = 1:numlocs, s = 1:numeffshifts, f = 1:numfragments[l,s]] >= 0, Int)	
    elseif lprelax_flag == 1
        @variable(ip, 0 <= x[i in orders, a in orderarcs.A[i]] <= 1)
	    @variable(ip, y[hasdriverarcs.A] >= 0)
	    @variable(ip, z[l = 1:numlocs, s = 1:numeffshifts, f = 1:numfragments[l,s]] >= 0)
    end
    @variable(ip, w[a in primaryarcs.A_space] >= 0)
	@variable(ip, ordtime[orders])
    @variable(ip, orderdelay[orders] >= 0)    #only used when deadlines turned on
  
	#Objective
	if deadlines_flag == 0
        @objective(ip, Min, lambda * sum((ordtime[i] - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in orders) + sum(c[a]*(y[a] ) for a in hasdriverarcs.A) + sum(u[a]*(w[a] ) for a in primaryarcs.A_space) )
    elseif deadlines_flag == 1
        @objective(ip, Min, lambda * sum(orderdelay[i] for i in orders) + sum(sum(c[a]*x[i,a] for a in orderarcs.A[i]) for i in orders) + sum(c[a]*(y[a] ) for a in hasdriverarcs.A) + sum(u[a]*(w[a] ) for a in primaryarcs.A_space) )
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
	#@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n],dummyarc)) for i in orders) - sum(y[a] for a in hasdriverarcs.A_plus[n]) == 0)
	@constraint(ip, truckFlowBalance[n in N_flow_t], sum(sum(x[i,a] for a in setdiff(orderarcs.A_minus[i,n], dummyarc)) for i in orders) + sum(y[a] for a in hasdriverarcs.A_minus[n]) - sum(sum(x[i,a] for a in setdiff(orderarcs.A_plus[i,n],dummyarc)) for i in orders) - sum(y[a] for a in hasdriverarcs.A_plus[n]) == 0)

	#Linking constraints
	@constraint(ip, driverAvailability[a in primaryarcs.A_space], sum(sum(sum(z[l,s,f] for f in fragmentscontaining[l,s,a]) for s in 1:numshifts) for l in 1:numlocs) == w[a]  )
	for i in orders, a in orderarcs.A_space[i]
		set_normalized_coefficient(driverAvailability[a], x[i,a], -1)
	end
	for a in hasdriverarcs.A_space
		set_normalized_coefficient(driverAvailability[a], y[a], -1)
	end

	#Driver constraints
	@constraint(ip, driverStartingLocs[l in 1:numlocs, s in 1:numeffshifts], sum(sum(z[l,s,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == length(driversets[l,s]))
	@constraint(ip, driverFlowBalance[l in 1:numlocs, s in 1:numeffshifts, n in N_flow_ls[l,s]], sum(z[l,s,f] for f in F_minus_ls[l,s,n]) - sum(z[l,s,f] for f in F_plus_ls[l,s,n]) == 0)

	optimize!(ip)

	ip_obj = objective_value(ip)
	if lprelax_flag == 1
		println("LP objective = ", ip_obj)
	else
		println("IP objective = ", ip_obj)
	end
	println("Time = ", solve_time(ip))

	totalmiles = sum(sum(c[a]*value(x[i,a]) for a in orderarcs.A[i]) for i in orders) + sum(c[a]*(value(y[a])) for a in hasdriverarcs.A) + sum(u[a]*(value(w[a]) ) for a in primaryarcs.A_space) 
	totaldelay = sum((value(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders)

	#Find the basis arcs
	orderArcSet_basis, orderArcSet_space_basis, A_plus_i_basis, A_minus_i_basis = getnonzeroarcs(value.(x), orderarcs)
	basisarcs = (A=orderArcSet_basis, A_space=orderArcSet_space_basis, A_minus=A_minus_i_basis, A_plus=A_plus_i_basis, available=[], closelocs=[]);
    
	return ip_obj, value.(x), value.(z), solve_time(ip), objective_bound(ip), basisarcs #, totalmiles, totaldelay
	
end

