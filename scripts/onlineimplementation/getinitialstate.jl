
function initializeonlinereporting()

	totalpastcost, totalordertrips, totalordermiles, ordermilesoutcomes, orderdelayoutcomes, totaldriverhours = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
	totalemptytrips, totalemptymiles, totaltaxitrips, totaltaximiles, ordermilespenalty, orderdelaypenalty = Dict(), Dict(), Dict(), Dict(), Dict(), Dict()
	total_delivtime, max_delivtime, shortestpossible_delivtime, shortestpossible_ordermiles = Dict(), Dict(), Dict(), Dict()
	for wklydelta in 0:timedelta:onlinetimehorizon+timedelta
		totalpastcost[wklydelta] = 0
		totalordertrips[wklydelta] = 0
		totalordermiles[wklydelta] = 0
		totaldriverhours[wklydelta] = 0
		totalemptytrips[wklydelta] = 0
		totalemptymiles[wklydelta] = 0
		totaltaxitrips[wklydelta] = 0
		totaltaximiles[wklydelta] = 0
		total_delivtime[wklydelta] = 0
		max_delivtime[wklydelta] = 0
		shortestpossible_delivtime[wklydelta] = 0
		shortestpossible_ordermiles[wklydelta] = 0
	end
	for i in currstate.orders
		ordermilesoutcomes[i] = 0
		orderdelayoutcomes[i] = 0
		ordermilespenalty[i] = 0
		orderdelaypenalty[i] = 0
	end
	pastordersegments, pastdriversegments_space, pastdriversegments, pastemptysegments, pasttaxisegments = [], [], [], [], []

	return total_delivtime, max_delivtime, shortestpossible_delivtime, shortestpossible_ordermiles, totalemptytrips, totalemptymiles, totaltaxitrips, totaltaximiles, ordermilespenalty, orderdelaypenalty, totalpastcost, totalordertrips, totalordermiles, ordermilesoutcomes, orderdelayoutcomes, totaldriverhours, pastordersegments, pastdriversegments_space, pastdriversegments, pastemptysegments, pasttaxisegments

end

#---------------------------------------------------------------------------------------#

function getinitialstate(nodesLookup, arcLookup, A_minus, A_plus, c, u)

    #Initialize orders
    orderwindowstart, orderwindowend = weekstart, weekstart + Dates.Hour(becomesavailablehours) - Dates.Second(1)
    orderwindowend = min(weekstart + Dates.Hour(horizon) + Dates.Millisecond(10), orderwindowend)
    includeorderidlist = generateorderlist(lhdataisbfilename, vntdataisbfilename, iterationordercap, numlocs)
    numorders, originloc, destloc, available, duedate, usedorderidlist, psseq, orderOriginalStartLoc, ordersinprogress = pullorders_initrivigoroutes(lhdataisbfilename, vntdataisbfilename, 10000, orderwindowstart, orderwindowend, tstep, horizon, prearcs, numlocs, tstepforordercreation, includeorderidlist)
    orders = [i for i in 1:numorders]
    highestorderindex = numorders
    Origin, Destination = formatorders(numorders, originloc, destloc, available, duedate, tstep)
    N_flow_i = flowbalancenodesets_i(orders, numnodes, Origin, Destination)
    orderOriginalStartTime, orderintransit_flag = findorderstartsandtransits(orders, Origin)

    #Initialize trucks
    m_0, m_end = truckdistribution(numtrucks, numlocs, N_0, N_end)
    loctruckcounter, trucksintransit = findtrucksintransit(ordersinprogress, originloc, available)
    m_0 = adjust_m_0(m_0, loctruckcounter)
    
    #Initialize drivers
    driversintransit, drivers, driverStartNodes, driverEndNodes, driverHomeLocs, assignedDrivers, N_flow_t, N_flow_d, alltimeswithinview, T_off_Monday8am, T_off, drivershift, T_off_0, T_off_constr, numshifts, T_on_0 = getdriverandshiftinfo()
    lasttimehome = zeros(length(drivers))

    #Distances
    distbetweenlocs, shortesttriptimes, shortestpatharclists, traveltimebetweenlocs_rdd, traveltimebetweenlocs_raw, traveltimebetweenlocs_llr = findtraveltimesanddistances(orders, Origin, Destination)
    nodesLookup, arcLookup, A_minus, A_plus, c, extendednodes, extendednumnodes, extendedarcs, extendednumarcs, u = extendtimespacenetwork(nodesLookup, arcLookup, A_minus, A_plus, c, u, distbetweenlocs)
    Destination = extendDestination(orders, Destination, extendednodes)

    #Define full state
    currstate = (m_0=m_0, m_end=m_end, trucksintransit=trucksintransit, orders=orders, 
    Origin=Origin, Destination=Destination, driversintransit=driversintransit, N_flow_d=N_flow_d, N_flow_i=N_flow_i,
    alltimeswithinview=alltimeswithinview, T_off=T_off, T_off_0=T_off_0, T_off_constr=T_off_constr, T_on_0=T_on_0,
    available=available, duedate=duedate, usedorderidlist=usedorderidlist, psseq=psseq, ordersinprogress=ordersinprogress, 
    shortesttriptimes=shortesttriptimes, orderintransit_flag=orderintransit_flag,
    driverStartNodes=driverStartNodes, driverEndNodes=driverEndNodes, assignedDrivers=assignedDrivers,
    lasttimehome=lasttimehome)

    return currstate, includeorderidlist, drivers, driverHomeLocs, drivershift, N_flow_t, T_off_Monday8am, numshifts, originloc, destloc, orderOriginalStartLoc, orderOriginalStartTime, highestorderindex, distbetweenlocs, shortestpatharclists, traveltimebetweenlocs_rdd, traveltimebetweenlocs_raw, traveltimebetweenlocs_llr, nodesLookup, arcLookup, A_minus, A_plus, c, u, extendednodes, extendednumnodes, extendedarcs, extendednumarcs

end

#---------------------------------------------------------------------------------------#