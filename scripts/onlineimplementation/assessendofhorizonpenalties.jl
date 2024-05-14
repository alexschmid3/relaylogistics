
function assessendofhorizonpenalties(currstate, currtime)

	for i in currstate.orders
		if orderOriginalStartTime[i] < onlinetimehorizon
            
            l1, l2 = nodesLookup[currstate.Origin[i][1]][1], nodesLookup[currstate.Destination[i][1]][1]

            #Add distance penalties
            totalpastcost[onlinetimehorizon+timedelta] += distbetweenlocs[l1,l2] * (1 + finallegdistancepenalty)
            totalordermiles[onlinetimehorizon+timedelta] += distbetweenlocs[l1,l2] * (1 + finallegdistancepenalty)
            shortestpossible_ordermiles[onlinetimehorizon+timedelta] += distbetweenlocs[originloc[i], destloc[i]]

            #Add distance penalty to outcomes for each order
            ordermilesoutcomes[i] += distbetweenlocs[l1,l2] * (1 + finallegdistancepenalty)
            ordermilespenalty[i] += distbetweenlocs[l1,l2] * (1 + finallegdistancepenalty)

            #Calculate delay penalty 
            posthorizonstarttime = max(0, orderOriginalStartTime[i] - currtime)
            if traveltimefordelay_flag == 0
                ttpenalty = traveltimebetweenlocs_rdd[l1,l2]
            elseif traveltimefordelay_flag >= 1
                ttpenalty = traveltimebetweenlocs_llr[l1,l2]
            end

            #Add delay penalties
            if orderOriginalStartTime[i] <= currtime
                orderdelay = ((currtime + ttpenalty * (1 + finallegtimepenalty) - orderOriginalStartTime[i]) - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i]
                orddeliverytime = currtime + ttpenalty * (1 + finallegtimepenalty) - orderOriginalStartTime[i]
            else
                orderdelay = ((ttpenalty * (1 + finallegtimepenalty)) - currstate.shortesttriptimes[i])/currstate.shortesttriptimes[i]
                orddeliverytime = ttpenalty * (1 + finallegtimepenalty) 
            end
            totalpastcost[onlinetimehorizon+timedelta] += lambda * orddeliverytime
            total_delivtime[onlinetimehorizon+timedelta] += orddeliverytime
            max_delivtime[onlinetimehorizon+timedelta] = max(max_delivtime[onlinetimehorizon+timedelta], orddeliverytime)
            shortestpossible_delivtime[onlinetimehorizon+timedelta] += currstate.shortesttriptimes[i]

            #Add segment to the list
			push!(pastordersegments, (i, onlinetimehorizon, dummyendtime, l1, l2))

            #Add delay penalty to outcomes for each order
            orderdelayoutcomes[i] += orddeliverytime
            orderdelaypenalty[i] += orddeliverytime

        end

	end

end
