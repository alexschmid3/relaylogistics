

function getnextorders(timedelta, currentdatetime, lh_data_file, vnt_data_file)

	#lh_data_file, vnt_data_file = "data/lh_data_isb_connect_clean.csv", "data/vnt_data_isb_connect_clean.csv"

	totaltimedelta = Dates.value(Dates.Hour(currentdatetime - weekstart))

	#====================================================#

	#Find new orders in the time horizon
	numexistingorders = length(currstate.orders)
	maxneworders = 10000 #iterationordercap #popfirst!(ordercaps)

	if maxneworders > 0

		#orderwindowstart, orderwindowend = weekstart + Dates.Hour(totaltimedelta) + Dates.Hour(becomesavailablehours) - Dates.Hour(timedelta) + Dates.Millisecond(10), weekstart + Dates.Hour(totaltimedelta) + Dates.Hour(becomesavailablehours)
		orderwindowstart = weekstart + Dates.Hour(totaltimedelta) + Dates.Hour(becomesavailablehours) - Dates.Hour(timedelta) + Dates.Millisecond(10)
		orderwindowend = min(weekstart + Dates.Hour(totaltimedelta) + Dates.Hour(becomesavailablehours), weekstart + Dates.Hour(onlinetimehorizon - 48) + Dates.Hour(becomesavailablehours))
		numneworders, originloc_ol, destloc_ol, available_ol, duedate_ol, orderidlist_new, psseq_ol, trueorderorigins = pullorders_rivigoroutes(lh_data_file, vnt_data_file, maxneworders, orderwindowstart, orderwindowend, currentdatetime + Dates.Hour(timedelta), tstep, horizon, prearcs, numlocs, timedelta, includeorderidlist)
		#newordersbyiter[totaltimedelta] = numneworders

		for orderid in orderidlist_new
			push!(currstate.usedorderidlist, orderid)
		end

		#Add to order list
        neworders = []
		for i in highestorderindex+1:highestorderindex+numneworders
			push!(currstate.orders, i)
            push!(neworders, i)
		end
        println(neworders)
		global highestorderindex += numneworders

		#Format Origin and Destination for each order
		currentindex = length(currstate.Origin)
		for i in 1:numneworders
			orderindex = currentindex + i
			currstate.Origin[orderindex] = []
			currstate.Destination[orderindex] = []
			for t in available_ol[i][1]:tstep:available_ol[i][2]
				push!(currstate.Origin[orderindex], nodes[originloc_ol[i], t])
			end
			for t in duedate_ol[i][1]:tstep:duedate_ol[i][2]
				push!(currstate.Destination[orderindex], nodes[destloc_ol[i], t])
			end
			push!(originloc, originloc_ol[i])
			push!(destloc, destloc_ol[i])
		end

		#Add pitstop sequences to psseq list
		for i in 1:numneworders
			orderindex = currentindex + i
			push!(currstate.psseq, psseq_ol[i])
		end

		#Update order flow balance node sets
		for i in 1:numneworders
			orderindex = currentindex + i
			currstate.N_flow_i[orderindex] = []
			for n in 1:numnodes
				if !(n in currstate.Origin[orderindex]) & !(n in currstate.Destination[orderindex])
					push!(currstate.N_flow_i[orderindex], n)
				end
			end
		end

		#Add order start times in terms of the original time horizon
		for i in 1:numneworders
			orderindex = currentindex + i
			orderOriginalStartTime[orderindex] = nodesLookup[currstate.Origin[orderindex][1]][2] + totaltimedelta
			push!(orderOriginalStartLoc, trueorderorigins[i] )
		end

		#Record shortest possible trip time for each order to be used in objective function
		for i in 1:numneworders
			orderindex = currentindex + i
			
			if traveltimefordelay_flag == 0
				shortestpathtime = traveltimebetweenlocs_rdd[nodesLookup[currstate.Origin[orderindex][1]][1], nodesLookup[currstate.Destination[orderindex][1]][1]]
			elseif traveltimefordelay_flag == 1
				shortestpathtime = traveltimebetweenlocs_raw[nodesLookup[currstate.Origin[orderindex][1]][1], nodesLookup[currstate.Destination[orderindex][1]][1]]
			elseif traveltimefordelay_flag == 2
				shortestpathtime = traveltimebetweenlocs_llr[nodesLookup[currstate.Origin[orderindex][1]][1], nodesLookup[currstate.Destination[orderindex][1]][1]]
			end

			push!(currstate.shortesttriptimes, shortestpathtime)
		end

		#====================================================#

        #Add dummy end destination for new orders
		for i in 1:numneworders
			orderindex = currentindex + i
			destinationlocation = nodesLookup[currstate.Destination[orderindex][1]][1]
			push!(currstate.Destination[orderindex], extendednodes[destinationlocation, dummyendtime])
		end

        #====================================================#
     
        #Find new order arcs and MAG arcs
        if operations == "relay"
            @time orderarcset_full, orderarcset_space_full, A_plus_i_full, A_minus_i_full = orderarcreduction(neworders, currstate.Origin, currstate.Destination)
        elseif operations == "ptp"
            @time orderarcset_full, orderarcset_space_full, A_plus_i_full, A_minus_i_full = pointotpointarcs(neworders, currstate.Origin, currstate.Destination)
        end
        if operations == "relay"
            @time magarcs = initializeorderarcsets(k, neworders, originloc, destloc, currstate.Origin, currstate.Destination, currstate.shortesttriptimes)
        elseif operations == "ptp"
            @time magarcs = (A=orderarcset_full, A_space=orderarcset_space_full, A_minus=A_minus_i_full, A_plus=A_plus_i_full, available=[])
        end

        #Update order arcs
        for i in neworders
            currarcs.orderarcs.A[i] = orderarcset_full[i]
            currarcs.orderarcs.A_space[i] = orderarcset_space_full[i]
            for n in 1:extendednumnodes
                currarcs.orderarcs.A_minus[i,n] = A_minus_i_full[i,n]
                currarcs.orderarcs.A_plus[i,n] = A_plus_i_full[i,n]
            end
        end

        #Update mag arcs
        for i in neworders
            currarcs.magarcs.A[i] = magarcs.A[i]
            currarcs.magarcs.A_space[i] = magarcs.A_space[i]
            for n in 1:extendednumnodes
                currarcs.magarcs.A_minus[i,n] = magarcs.A_minus[i,n]
                currarcs.magarcs.A_plus[i,n] = magarcs.A_plus[i,n]
            end
        end

        for i in neworders
            ordermilesoutcomes[i] = 0
            ordermilespenalty[i] = 0
            orderdelayoutcomes[i] = 0
            orderdelaypenalty[i] = 0
        end

	else

        1+1
		#newordersbyiter[totaltimedelta] = 0

	end

end
