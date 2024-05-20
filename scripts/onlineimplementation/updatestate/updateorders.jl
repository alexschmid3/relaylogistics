
function updateorders(x, timedelta, currentdatetime, basisarcs)

	#Update the time 
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

	#Find current location of orders from last time (may be in transit) and update origin/destination windows
	completedorders = []
    orderintransit_flag = Dict()
	for i in currstate.orders
		orderintransit_flag[i] = 0 
	end
	
	for i in currstate.orders
		originstarttime = nodesLookup[currstate.Origin[i][1]][2]

		#If order path infeasible in current iteration
		if value(x[i, dummyarc]) > 0.001
			#Iterate all origin nodes forward by timedelta
			originloc2, originaltime = nodesLookup[Origin[i][1]]
			newstarttime = max(0,  originaltime - timedelta)

			#If order in progress, move the origin node forward by timedelta
			if i in ordersinprogress
				newnode = nodes[originloc2, newstarttime]
				Origin[i] = [newnode]

			#If order not in progress, move start of origin forward by timedelta, but leaving end of origin as the horizon
			else
				for t in originaltime-tstep:-tstep:newstarttime
					newnode = nodes[originloc2, t]
					pushfirst!(Origin[i], newnode)
				end
			end
		
			#Iterate all destination nodes forward by timedelta
			newDestinationList = []
			n = currstate.Destination[i][1]
			destlocation, originaltime = nodesLookup[n]
			for t in max(0,originaltime-timedelta):tstep:horizon
				newnode = nodes[destlocation, t]
				push!(newDestinationList, newnode)
			end
			push!(newDestinationList, extendednodes[destlocation, dummyendtime])
			currstate.Destination[i] = newDestinationList

		#Else if the order was available before the start of the new time horizon
		elseif originstarttime < timedelta

			#Find current location of order
			#Initialize with origin location and time
			currentloc, availtime, endwindowtime, mostrecentarc = nodesLookup[currstate.Origin[i][1]][1], 0.0, nodesLookup[last(currstate.Origin[i])][2], 0
			
            #println("$i, $availtime, $endwindowtime")

            #Find the arcs used by the order in the current solution
            usedarcs = [a for a in orderarcs.A[i] if value(x[i,a]) > 1e-4]
            usedarcs = sort(usedarcs, by=x->nodesLookup[arcLookup[x][1]][2])

            #Find the updated order location and time
            mostrecentarc = 0
            for a in usedarcs
                l1,t1 = nodesLookup[arcLookup[a][1]]
                l2,t2 = nodesLookup[arcLookup[a][2]]
                if t1 < timedelta - 1e-4
                    currentloc, availtime = l2,t2 
                    #println("($i,$a), $availtime, $currentloc")
                    mostrecentarc = a
                else 
                    break
                end
            end
            availtime = availtime == 0 ? timedelta : availtime
            orderintransit_flag[i] = availtime > timedelta ? 1 : 0
            endwindowtime = currentloc == orderOriginalStartLoc[i] ? horizon+timedelta : availtime

            #println("$i, $availtime, $endwindowtime")

            #Move forward by timedelta
            availtime = availtime - timedelta
            endwindowtime = endwindowtime - timedelta

            #println("$i, $availtime, $endwindowtime")

			#Check whether order was completed
			if currentloc == nodesLookup[currstate.Destination[i][1]][1]

                push!(completedorders, i)
                remove!(currstate.ordersinprogress, i)

				arcendtime = arcfinishtime[mostrecentarc] - timedelta
                orddeliverytime = arcendtime + totaldelta - orderOriginalStartTime[i]
				orddelay = ((arcendtime + totaldelta - orderOriginalStartTime[i]) - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i]
				totalpastcost[totaldelta-timedelta] += lambda * orddeliverytime
				total_delivtime[totaldelta-timedelta] += orddeliverytime
				max_delivtime[totaldelta-timedelta] = max(max_delivtime[totaldelta-timedelta], orddeliverytime)		
				orderdelayoutcomes[i] += arcendtime + totaldelta - orderOriginalStartTime[i]
                shortestpossible_ordermiles[totaldelta-timedelta] += distbetweenlocs[originloc[i], destloc[i]]
                shortestpossible_delivtime[totaldelta-timedelta] += currstate.shortesttriptimes[i]
			else
				#Create new origin set
				newOriginList = []
				for t in availtime:tstep:endwindowtime
					newnode = nodes[currentloc, t]
					push!(newOriginList, newnode)
				end
				currstate.Origin[i] = newOriginList

				#Mark if the order has left its destination
				if currentloc != orderOriginalStartLoc[i]
					#Add to ordersinprogress
					if !(i in currstate.ordersinprogress)
						push!(currstate.ordersinprogress, i)
					end
				end

				#Iterate all destination nodes forward by timedelta
				newDestinationList = []
				n = currstate.Destination[i][1]
				destlocation, originaltime = nodesLookup[n]
				for t in max(0,originaltime-timedelta):tstep:horizon
					newnode = nodes[destlocation, t]
					push!(newDestinationList, newnode)
				end
				push!(newDestinationList, extendednodes[destlocation, dummyendtime])
				currstate.Destination[i] = newDestinationList
				
			end
			
		#Else, the order becomes available after the start of the new time horizon
		else
			#Iterate all origin nodes forward by timedelta
			originloc2, originaltime = nodesLookup[currstate.Origin[i][1]]
			newstarttime = max(0,  originaltime - timedelta)

			#If order in progress, move the origin node forward by timedelta
			if i in currstate.ordersinprogress
				newnode = nodes[originloc2, newstarttime]
				currstate.Origin[i] = [newnode]

			#If order not in progress, move start of origin forward by timedelta, but leaving end of origin as the horizon
			else
				for t in originaltime-tstep:-tstep:newstarttime
					newnode = nodes[originloc2, t]
					pushfirst!(currstate.Origin[i], newnode)
				end
			end

			#Iterate all destination nodes forward by timedelta
			newDestinationList = []
			n = currstate.Destination[i][1]
			destlocation, originaltime = nodesLookup[n]
			for t in max(0,originaltime-timedelta):tstep:horizon
				newnode = nodes[destlocation, t]
				push!(newDestinationList, newnode)
			end
			push!(newDestinationList, extendednodes[destlocation, dummyendtime])
			currstate.Destination[i] = newDestinationList
		end
	end

	#====================================================#

	#Remove completed orders
    println("Completed ", length(completedorders), " orders")
	for i in completedorders
		remove!(currstate.orders, i)
	end
 
	#====================================================#

	#Update order flow balance node sets
	for i in currstate.orders
		currstate.N_flow_i[i] = []
		for n in 1:numnodes
			if !(n in currstate.Origin[i]) & !(n in currstate.Destination[i])
				push!(currstate.N_flow_i[i], n)
			end
		end
	end

end
