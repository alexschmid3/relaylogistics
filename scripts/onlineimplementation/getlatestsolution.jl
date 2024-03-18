
function updatepastsegments(timedelta, x, y, z, w, restrictedArcSet, currfragments, driversingroup)

	wklydelta = mod(Dates.value(Dates.Hour(currentdatetime - weekstart)), 168)
	
	#====================================================#

	#Find the set of arcs that are locked in after the next iteration
	lockedarcs = []
	for a in 1:numarcs
		if nodesLookup[arcLookup[a][1]][2] < timedelta
			push!(lockedarcs, a)
		end
	end
	lockedjourneys = Dict()
	for (hl,ss,sn,aln) in currfragments.driversets
		lockedjourneys[hl,ss,sn,aln] = []
		for l in 1:numlocs, t in 0:tstep:timedelta-tstep
			lockedjourneys[hl,ss,sn,aln] = union(lockedjourneys[hl,ss,sn,aln], F_plus_g[hl,ss,sn,aln,nodes[l,t]])
		end
	end

	#====================================================#

	#Add order segments
	for i in orders, a in intersect(lockedarcs, setdiff(restrictedArcSet[i], dummyarc))
		if getvalue(x[i,a]) > 0.001
			#Add new current segment = (i, starttime, endtime, startloc, endloc)
			arcstartloc, arcendloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
			push!(pastordersegments, (i, nodesLookup[arcLookup[a][1]][2] + wklydelta, nodesLookup[arcLookup[a][2]][2] + wklydelta, arcstartloc, arcendloc))
			#Update the total past cost with the completed segment
			global totalpastcost += c[a]
			global totalordertrips += 1
			global totalordermiles += c[a]
			ordermilesoutcomes[i] += c[a]
		end
	end

	#====================================================#

	#Identify and connect journeys
	driverarcstaken = Dict()
	for d in drivers
		driverarcstaken[d] = []
	end
	for (hl,ss,sn,aln) in driversets
		driverstoassign = copy(driversingroup[hl,ss,sn,aln])
		usedjourneys = []
		for f in lockedjourneys[hl,ss,sn,aln]
			numdriversonfragment = convert(Int64, round(getvalue(z[(hl,ss,sn,aln,f)]),digits=0))
			for fragcopy in 1:numdriversonfragment
				push!(usedjourneys, f)
			end
		end

		while driverstoassign != []
			d = pop!(driverstoassign)	
			currnode = sn
			currfrag = -1
			while [f for f in usedjourneys if f in F_plus_g[hl,ss,sn,aln,currnode]] != []
				#Choose next fragment and remove from the available fragment list
				currfrag = pop!([f for f in usedjourneys if f in F_plus_g[hl,ss,sn,aln,currnode]])
				currfragindex = findfirst(x->x==currfrag, usedjourneys)
				deleteat!(usedjourneys, currfragindex)
				
				#Add arcs on the fragment
				addarclist = sort(intersect(fragmentarcs[hl,ss,sn,aln,currfrag], lockedarcs), by=x->nodesLookup[arcLookup[x][2]][2])
				for a in addarclist
					arcstartloc, arcendloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
					push!(driverarcstaken[d], a)
					push!(pastdriversegments, (d, nodesLookup[arcLookup[a][1]][2] + wklydelta, nodesLookup[arcLookup[a][2]][2] + wklydelta, arcstartloc, arcendloc))
					if a in A_space
						push!(pastdriversegments_space, (d, nodesLookup[arcLookup[a][1]][2] + wklydelta, nodesLookup[arcLookup[a][2]][2] + wklydelta, nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
					end
					if arcstartloc != arcendloc
						global totaldriverhours += arcLength_raw[arcstartloc, arcendloc]
					end
				end

				#Next node
				lastarc = last(addarclist)
				currnode = arcLookup[lastarc][2]
			end
		end
	end

    #Add empty segments
	for a in intersect(lockedarcs, A_hasdriver_space)
		if getvalue(y[a]) > 0.001
			#Add new current segment = (numtrucks, starttime, endtime, startloc, endloc)
			push!(pastemptysegments, (getvalue(y[a]), nodesLookup[arcLookup[a][1]][2] + wklydelta, nodesLookup[arcLookup[a][2]][2] + wklydelta, nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
			global totalpastcost += c[a] * getvalue(y[a])
			global totalemptytrips += getvalue(y[a])
			global totalemptymiles += c[a] * getvalue(y[a])
		end
	end

	#Add taxi segments
	for a in intersect(lockedarcs, A_space)
		if getvalue(w[a]) > 0.001
			#Add new current segment = (numtrucks, starttime, endtime, startloc, endloc)
			push!(pasttaxisegments, (getvalue(w[a]), nodesLookup[arcLookup[a][1]][2] + wklydelta, nodesLookup[arcLookup[a][2]][2] + wklydelta, nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
			global totalpastcost += u[a] * getvalue(w[a])
			global totaltaxitrips += getvalue(w[a])
			global totaltaximiles += u[a] * getvalue(w[a])
		end
	end

	return driverarcstaken

end