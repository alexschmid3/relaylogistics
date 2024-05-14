
function updatetrucks(timedelta, currentdatetime, weekstart, x, y)

	#Update the time 
	wklydelta = mod(Dates.value(Dates.Hour(currentdatetime - weekstart)), 168)

	#====================================================#

	#Create set of arcs that begin before the time horizon
	earlyarcs = []
	for a in 1:numarcs
		if (nodesLookup[arcLookup[a][1]][2] < timedelta) & (nodesLookup[arcLookup[a][2]][2] > timedelta)
			push!(earlyarcs, a)
		end
	end

	#Find any trucks that are currently in transit
	intransittruckcount = [0 for l in 1:numlocs]
    for item in currstate.trucksintransit
        remove!(currstate.trucksintransit, item)
    end
	for t in timedelta:tstep:horizon, l in 1:numlocs
		n = nodes[l, t]
		totaltrucks = 0
		for i in currstate.orders
			if intersect(earlyarcs, currarcs.magarcs.A_minus[i,n]) != []
				totaltrucks += sum(value(x[i,a]) for a in intersect(earlyarcs, currarcs.magarcs.A_minus[i,n]))
			end
		end
		if intersect(earlyarcs, currarcs.hasdriverarcs.A_minus[n]) != []
			totaltrucks += sum(value(y[a]) for a in intersect(earlyarcs, currarcs.hasdriverarcs.A_minus[n])) 
		end
		totaltrucks = round(totaltrucks)
		#totaltrucks = sum(sum(value(x[i,a]) for a in intersect(earlyarcs, currarcs.magarcs.A_minus[i,n])) for i in orders) + sum(value(y[a]) for a in intersect(earlyarcs, A_minus_hd[n])) 
		if totaltrucks > 0.001
			intransittruckcount[l] += totaltrucks
			push!(currstate.trucksintransit, (l, t - timedelta, totaltrucks))
		end
	end

	#Update truck initial locations
	for l in 1:numlocs
		n = nodes[l, timedelta]
		enteringtrucks = 0 
		for i in currstate.orders
			if setdiff(currarcs.magarcs.A_minus[i,n], dummyarc) != []
				enteringtrucks += sum(value(x[i,a]) for a in setdiff(currarcs.magarcs.A_minus[i,n], dummyarc)) 
			end
		end
		if currarcs.hasdriverarcs.A_minus[n] != []
			enteringtrucks += sum(value(y[a]) for a in currarcs.hasdriverarcs.A_minus[n]) 
		end
		#enteringtrucks += sum(sum(value(x[i,a]) for a in setdiff(A_minus_i[i,n], dummyarc)) for i in orders) + sum(value(y[a]) for a in A_minus_hd[n])
		enteringtrucks = round(enteringtrucks)
		#Total trucks starting at location l = the entering trucks plus any trucks found to be in transit to location l
		currstate.m_0[l] = enteringtrucks + intransittruckcount[l]
	end

	println("Total m_0 = ", sum(values(currstate.m_0)) )
	@assert(sum(values(currstate.m_0)) == numtrucks)

end

