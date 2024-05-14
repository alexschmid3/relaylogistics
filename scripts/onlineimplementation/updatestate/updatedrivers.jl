
#-------------------------------------------------------------------------------------#

function updatelasttimehome(driverarcstaken)

	for d in drivers
        #Update time based on arcs
        for a in sort(driverarcstaken[d], by=x->nodesLookup[arcLookup[x][1]][2])
            arcendloc, arcendtime = nodesLookup[arcLookup[a][2]][1], nodesLookup[arcLookup[a][2]][2]
			if arcendloc == driverHomeLocs[d]
                currstate.lasttimehome[d] = arcendtime
            end
        end

        #Move forward in time
        currstate.lasttimehome[d] = currstate.lasttimehome[d] - timedelta
	end

end

#-------------------------------------------------------------------------------------#

function updatedriverlocations(currentdatetime, driverarcstaken)

	wklydelta = mod(Dates.value(Dates.Hour(currentdatetime - weekstart)), 168)

	#====================================================#

	#Update driver start nodes and find drivers in transit
	for item in currstate.driversintransit
        remove!(currstate.driversintransit, item)
    end
	for d in drivers
		if driverarcstaken[d] != []
            #println("$d --> ", driverarcstaken[d])
			driverarcstaken[d] = sort!(driverarcstaken[d], by=x->nodesLookup[arcLookup[x][2]][2])
			newstartloc, newstarttime = nodesLookup[arcLookup[last(driverarcstaken[d])][2]][1], nodesLookup[arcLookup[last(driverarcstaken[d])][2]][2] - timedelta
			currstate.driverStartNodes[d] = nodes[newstartloc, newstarttime]
			if newstarttime > 0
				push!(currstate.driversintransit, (d, newstartloc, newstarttime))
			end	
		else
			newstartloc, newstarttime = nodesLookup[currstate.driverStartNodes[d]][1], nodesLookup[currstate.driverStartNodes[d]][2] - timedelta
			currstate.driverStartNodes[d] = nodes[newstartloc, newstarttime]
			if (newstarttime > 0) & (newstarttime < horizon)
				push!(currstate.driversintransit, (d, newstartloc, newstarttime))
			end	
		end	
	end

	#====================================================#

	#Update driver flow balance node set
	for d in drivers
		currstate.N_flow_d[d] = setdiff([n for n in 1:numnodes], union(N_end, currstate.driverStartNodes[d]))
	end

	#====================================================#

	#Count nights away
	#=for d in drivers, a in driverarcstaken[d]
		l, t = nodesLookup[arcLookup[a][1]]
		if (currtime+t in T_off_Monday8am_0[drivershift[d]]) & (l == driverHomeLocs[d]) 
			drivernightshome[d] -= 1
			drivernightsaway[d] += 1
		elseif (solutionmethod == "ptp") & (a in A_space) #Drivers may be "moving" on a long arc during their off hours
			for t2 in t:tstep:nodesLookup[arcLookup[a][2]][2]
				if currtime+t2 in T_off_Monday8am_0[drivershift[d]]
					drivernightshome[d] -= 1
					drivernightsaway[d] += 1
				end
			end
		end
	end=#

end

#-------------------------------------------------------------------------------------#

function updatedriversshifts(currentdatetime, weekstart, T_off_Monday8am)

	#Update the time 
	wklydelta = mod(Dates.value(Dates.Hour(currentdatetime - weekstart)), 168)

	#====================================================#

	#UPDATE DRIVER SCHEDULES

	#Adjust driver shift schedule for the new date and time
	for ss in 1:length(T_off_Monday8am)
		shiftedshift = []
		for hr in T_off_Monday8am[ss]
			if wklydelta <= hr <= wklydelta+horizon + 2*24*maxnightsaway
				push!(shiftedshift, hr - wklydelta)
			end
		end
		currstate.T_off[ss] = shiftedshift
	end

	#Create other driver shift lists (used to write constraints)
	for d in drivers
		currstate.T_off_0[d], currstate.T_off_constr[d] = [], []
		for t in 0:tstep:horizon + 2*24*maxnightsaway
			if (t in currstate.T_off[drivershift[d]]) & !(t-tstep in currstate.T_off[drivershift[d]])
				push!(currstate.T_off_0[d], t)
				if (t + 24 <= horizon) & (intersect([t3 for t3 in t:tstep:t+24],currstate.T_off_0[d]) != [])
					push!(currstate.T_off_constr[d], t)  
				end
			end
		end
	end
	for d in drivers
		currstate.T_on_0[d] = []
		for t in 0:tstep:horizon+24*2*maxnightsaway
			if (t in setdiff(0:tstep:horizon+2*24*maxnightsaway, currstate.T_off[drivershift[d]])) & !(t-tstep in setdiff(0:tstep:horizon+2*24*maxnightsaway, currstate.T_off[drivershift[d]]))
				push!(currstate.T_on_0[d], t)
			end
		end
	end

end

#-------------------------------------------------------------------------------------#