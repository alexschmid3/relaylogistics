

function updatepastsegments(timedelta, x, y, z, w, candidatejourneys, currentdatetime, basisarcs)

	# x, y, z, w = x_ip, y_ip, z_ip, w_ip

    wklydelta = mod(Dates.value(Dates.Hour(currentdatetime - weekstart)), 168)
    totaldelta = Dates.value(Dates.Hour(currentdatetime - weekstart))

    #====================================================#

    if (operations == "relay") & (solutionmethod == "mag")
		orderarcs = currarcs.magarcs
	elseif (operations == "relay") & (solutionmethod == "basisip")
		orderarcs = basisarcs
	else	
		orderarcs = currarcs.orderarcs
	end

    #====================================================#

    #Find the set of arcs that are locked in after the next iteration
    lockedarcs = []
    for a in 1:numarcs
        if nodesLookup[arcLookup[a][1]][2] < timedelta
            push!(lockedarcs, a)
        end
    end
    lockedjourneys = Dict()
    if candidatejourneys == -1
        for (hl,ss,sn,lth) in currfragments.driversets
            lockedjourneys[hl,ss,sn,lth] = []
            for l in 1:numlocs, t in 0:tstep:timedelta-tstep
                lockedjourneys[hl,ss,sn,lth] = union(lockedjourneys[hl,ss,sn,lth], currfragments.F_plus_g[hl,ss,sn,lth,nodes[l,t]])
            end
        end
    else
        for (hl,ss,sn,lth) in currfragments.driversets
            lockedjourneys[hl,ss,sn,lth] = []
            for l in 1:numlocs, t in 0:tstep:timedelta-tstep
                lockedjourneys[hl,ss,sn,lth] = union(lockedjourneys[hl,ss,sn,lth], intersect(candidatejourneys[hl,ss,sn,lth],currfragments.F_plus_g[hl,ss,sn,lth,nodes[l,t]]))
            end
        end
    end
    #====================================================#

    #Add order segments
    for i in currstate.orders, a in intersect(lockedarcs, setdiff(orderarcs.A[i], dummyarc))
        if value(x[i,a]) > 0.001
            #Add new current segment = (i, starttime, endtime, startloc, endloc)
            arcstartloc, arcendloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
            push!(pastordersegments, (i, nodesLookup[arcLookup[a][1]][2] + totaldelta, nodesLookup[arcLookup[a][2]][2] + totaldelta, arcstartloc, arcendloc))
            #Update the total past cost with the completed segment
            totalpastcost[totaldelta] += c[a]
            totalordertrips[totaldelta] += 1
            totalordermiles[totaldelta] += c[a]
            ordermilesoutcomes[i] += c[a]
        end
    end

    #====================================================#

    #Identify and connect journeys
    driverarcstaken = Dict()
    for d in drivers
        driverarcstaken[d] = []
    end
    for (hl,ss,sn,lth) in currfragments.driversets
        driverstoassign = copy(currfragments.driversingroup[hl,ss,sn,lth])
        usedjourneys = []
        for f in lockedjourneys[hl,ss,sn,lth]
            numdriversonfragment = convert(Int64, round(value(z[(hl,ss,sn,lth),f]),digits=0))
            for fragcopy in 1:numdriversonfragment
                push!(usedjourneys, f)
            end
        end
        
        while driverstoassign != []
            d = pop!(driverstoassign)	
            currnode = sn
            currfrag = -1
            while [f for f in usedjourneys if f in currfragments.F_plus_g[hl,ss,sn,lth,currnode]] != []
                #Choose next fragment and remove from the available fragment list
                currfrag = pop!([f for f in usedjourneys if f in currfragments.F_plus_g[hl,ss,sn,lth,currnode]])
                currfragindex = findfirst(x->x==currfrag, usedjourneys)
                deleteat!(usedjourneys, currfragindex)
                
                #Add arcs on the fragment
                addarclist = sort(intersect(currfragments.fragmentarcs[hl,ss,sn,lth,currfrag], lockedarcs), by=x->nodesLookup[arcLookup[x][2]][2])
                for a in addarclist
                    arcstartloc, arcendloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
                    push!(driverarcstaken[d], a)
                    push!(pastdriversegments, (d, nodesLookup[arcLookup[a][1]][2] + totaldelta, nodesLookup[arcLookup[a][2]][2] + totaldelta, arcstartloc, arcendloc))
                    if a in A_space
                        push!(pastdriversegments_space, (d, nodesLookup[arcLookup[a][1]][2] + totaldelta, nodesLookup[arcLookup[a][2]][2] + totaldelta, nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
                    end
                    if arcstartloc != arcendloc
                        totaldriverhours[totaldelta] += arcLength_raw[arcstartloc, arcendloc]
                    end
                end

                #Next node
                if addarclist == []
                    break
                else
                    lastarc = last(addarclist)
                    currnode = arcLookup[lastarc][2]
                end
            end
        end
    end

    #Add empty segments
    for a in intersect(lockedarcs, currarcs.hasdriverarcs.A, primaryarcs.A_space)
        if value(y[a]) > 0.001
            #Add new current segment = (numtrucks, starttime, endtime, startloc, endloc)
            push!(pastemptysegments, (value(y[a]), nodesLookup[arcLookup[a][1]][2] + totaldelta, nodesLookup[arcLookup[a][2]][2] + totaldelta, nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
            totalpastcost[totaldelta] += c[a] * value(y[a])
            totalemptytrips[totaldelta] += value(y[a])
            totalemptymiles[totaldelta] += c[a] * value(y[a])
        end
    end

    #A_space_all = primaryarcs.A_space
    #for i in currstate.orders
    #    goodones = [a for a in orderarcs.A_space[i] if nodesLookup[arcLookup[a][1]][2] < horizon]
    #    A_space_all = union(A_space_all, goodones)
    #end
    A_space_all = []
	for a in 1:extendednumarcs
		l1,t1 = nodesLookup[arcLookup[a][1]]
		l2,t2 = nodesLookup[arcLookup[a][2]]
		if (l1 != l2) & (t1 < horizon)
			push!(A_space_all, a)
		end
	end

    #Add taxi segments
    for a in intersect(lockedarcs, A_space_all)
        if value(w[a]) > 0.001
            #Add new current segment = (numtrucks, starttime, endtime, startloc, endloc)
            push!(pasttaxisegments, (value(w[a]), nodesLookup[arcLookup[a][1]][2] + totaldelta, nodesLookup[arcLookup[a][2]][2] + totaldelta, nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
            totalpastcost[totaldelta] += u[a] * value(w[a])
            totaltaxitrips[totaldelta] += value(w[a])
            totaltaximiles[totaldelta] += u[a] * value(w[a])
        end
    end

	return driverarcstaken

end

