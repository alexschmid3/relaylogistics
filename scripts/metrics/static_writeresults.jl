
#---------------------------------------------------------------------------------------#
#-------------------------------------WRITE RESULTS-------------------------------------#  
#---------------------------------------------------------------------------------------#

function writeresults_abcg_static(filename, cg_iter, rmp_obj, rmp_time, pp_time, pp_time_par, x, y, z, w, ordtime, currtime, minreducedcost)

	bigOrderArcSet = []
	for i in orders
       bigOrderArcSet = union(bigOrderArcSet, orderArcSet[i])
    end

	#Calculate necessary quantities
	totalorders = length(orders)

	obj_delay = sum((getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) 
	obj_dist = sum(sum(c[a]*getvalue(x[i,a]) for a in orderArcSet[i]) for i in orders) + sum(c[a]*getvalue(y[a]) for a in A_hasdriver) + sum(u[a]*getvalue(w[a]) for a in A_space) 
	infeasible_order_count = sum(getvalue(x[i,dummyarc]) for i in orders)
	if intersect(bigOrderArcSet, [a for a in numarcs+1:extendednumarcs]) == []
		incomplete_order_count = 0
	else
		incomplete_order_count = sum(sum(getvalue(x[i,a]) for a in intersect(orderArcSet[i], [a for a in numarcs+1:extendednumarcs])) for i in orders if intersect(orderArcSet[i], [a for a in numarcs+1:extendednumarcs]) != [])
	end

	delay_per_order = obj_delay/length(orders)
	max_order_delay = maximum([(getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders]) 
	totalnightsathome = sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)
	totalnightsaway = sum(length(T_off_0[d]) for d in drivers) - sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)

	totalordertrips, totalemptytrips, totaltaxitrips = 0, 0, 0
	totalordermiles, totalemptymiles, totaltaximiles = 0, 0, 0
	for i in orders, a in intersect(A_space, orderArcSet[i])
		if getvalue(x[i,a]) > 0.001
			totalordertrips += getvalue(x[i,a])
			totalordermiles += c[a] *getvalue(x[i,a])
		end
	end
	for a in intersect(A_space, A_hasdriver)
		if getvalue(y[a]) > 0.001
			totalemptytrips += getvalue(y[a])
			totalemptymiles += c[a]*getvalue(y[a])
		end
	end

	for a in A_space
		if getvalue(w[a]) > 0.001
			totaltaxitrips += getvalue(w[a])
			totaltaximiles += u[a]*getvalue(w[a])
		end
	end

	#----------------------------WRITE ORDER TIMES TO CSV----------------------------#

	df = DataFrame(ID = [runid],
			example = [ex],
			timehorizon = [horizon], 
			timestep = [tstep], 
			numlocs = [numlocs], 
			numdrivers = [length(drivers)], 
			numtrucks = [numtrucks], 
			numarcs = [numarcs], 
			lambdaobj = [lambda],
			numorders = [totalorders], 
			abcgiteration = [cg_iter],
			heuristicdp = [solvedpheuristically_flag_now],
			totalabcgiteration = [""],
			minreducedcost = [minreducedcost],
			optval = [rmp_obj], 
			optval_delay = [obj_delay], 
			optval_dist = [obj_dist], 
			infeasible_order_count = [infeasible_order_count],
			incomplete_order_count = [incomplete_order_count],
			ordertrips = [totalordertrips], 
			emptytrips = [totalemptytrips], 
			taxitrips = [totaltaxitrips],
			ordermiles = [totalordermiles], 
			emptymiles = [totalemptymiles], 
			taximiles = [totaltaximiles], 
			delay_per_order = [delay_per_order], 
			max_order_delay = [max_order_delay],
			nightsathome = [totalnightsathome],
			nightsaway = [totalnightsaway],
			ip_time = [""],
			rmp_time = [rmp_time],
			pp_time = [pp_time],
			pptime_parallel = [pp_time_par],
			order_arc_count = [sum(length(orderArcSet[i]) for i in orders)]
           )

	if cg_iter == 1
		CSV.write(filename, df)
	else
		CSV.write(filename, df, append=true)
	end

end

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function writeresults_ip_abcg_static(filename, comptimesfilename, total_cg_iter, ip_obj, ip_time, cg_rmptimes, cg_pptimes, cg_pptimes_par, x, y, z, w, ordtime)

	bigOrderArcSet = []
	for i in orders
       bigOrderArcSet = union(bigOrderArcSet, orderArcSet[i])
    end

	#Calculate necessary quantities
	totalorders = length(orders)

	obj_delay = sum((getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) 
	obj_dist = sum(sum(c[a]*getvalue(x[i,a]) for a in orderArcSet[i]) for i in orders) + sum(c[a]*getvalue(y[a]) for a in A_hasdriver) + sum(u[a]*getvalue(w[a]) for a in A_space) 
	infeasible_order_count = sum(getvalue(x[i,dummyarc]) for i in orders)
	if intersect(bigOrderArcSet, [a for a in numarcs+1:extendednumarcs]) == []
		incomplete_order_count = 0
	else
		incomplete_order_count = sum(sum(getvalue(x[i,a]) for a in intersect(orderArcSet[i], [a for a in numarcs+1:extendednumarcs])) for i in orders if intersect(orderArcSet[i], [a for a in numarcs+1:extendednumarcs]) != [])
	end

	delay_per_order = obj_delay/length(orders)
	max_order_delay = maximum([(getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders]) 
	totalnightsathome = sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)
	totalnightsaway = sum(length(T_off_0[d]) for d in drivers) - sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)

	totalordertrips, totalemptytrips, totaltaxitrips = 0, 0, 0
	totalordermiles, totalemptymiles, totaltaximiles = 0, 0, 0
	for i in orders, a in intersect(A_space, orderArcSet[i])
		if getvalue(x[i,a]) > 0.001
			totalordertrips += 1
			totalordermiles += c[a]
		end
	end
	for a in intersect(A_space, A_hasdriver)
		if getvalue(y[a]) > 0.001
			totalemptytrips += getvalue(y[a])
			totalemptymiles += c[a]*getvalue(y[a])
		end
	end

	for a in A_space
		if getvalue(w[a]) > 0.001
			totaltaxitrips += getvalue(w[a])
			totaltaximiles += u[a]*getvalue(w[a])
		end
	end

	rmp_time = sum(cg_rmptimes)
	cg_pp_time = sum(cg_pptimes)
	cg_pptime_parallel = sum(cg_pptimes_par)

	orderarccnt = sum(length(orderArcSet[i]) for i in orders)

	#----------------------------SAVE TIMES FOR REPORTING----------------------------#

	#push!(totalcgiterlist, total_cg_iter)
	#push!(iptimelistlist, ip_time)
	#push!(rmptimelist, rmp_time)
	#push!(cgpptimelist, cg_pp_time)
	#push!(orderarccountlist, orderarccnt)

	#----------------------------WRITE ORDER TIMES TO CSV----------------------------#

	df = DataFrame(ID = [runid],
			example = [ex],
			timehorizon = [horizon], 
			timestep = [tstep], 
			numlocs = [numlocs], 
			numdrivers = [length(drivers)], 
			numtrucks = [numtrucks], 
			numarcs = [numarcs], 
			lambdaobj = [lambda],
			numorders = [length(orders)], 
			abcgiteration = ["IP"],
			heuristicdp = [solvedpheuristically_flag_now],
			totalabcgiteration = [total_cg_iter],
			minreducedcost = [""],
			optval = [ip_obj], 
			optval_delay = [obj_delay], 
			optval_dist = [obj_dist], 
			infeasible_order_count = [infeasible_order_count],
			incomplete_order_count = [incomplete_order_count],
			ordertrips = [totalordertrips], 
			emptytrips = [totalemptytrips], 
			taxitrips = [totaltaxitrips],
			ordermiles = [totalordermiles], 
			emptymiles = [totalemptymiles], 
			taximiles = [totaltaximiles], 
			delay_per_order = [delay_per_order], 
			max_order_delay = [max_order_delay],
			nightsathome = [totalnightsathome],
			nightsaway = [totalnightsaway],
			ip_time = [ip_time],
			cg_rmp_time = [rmp_time],
			cg_pp_time = [cg_pp_time],
			cg_pptime_parallel = [cg_pptime_parallel],
			order_arc_count = [orderarccnt]
           )

	CSV.write(filename, df, append=true)

	#--------------------------------------------------------------------------------#

	df2 = DataFrame(ID = [runid],
			example = [ex],
			numorders = [length(orders)], 
			numlocs = [numlocs], 
			numdrivers = [length(drivers)], 
			numtrucks = [numtrucks], 
			timehorizon = [horizon], 
			timestep = [tstep], 
			lambdaobj = [lambda],
			objval = [ip_obj],
			delayobj = [obj_delay],
			milesobj = [obj_dist],
			totaltime = [ip_time + rmp_time + cg_pptime_parallel],
			order_arc_count = [orderarccnt]
           )

	CSV.write(comptimesfilename, df2)

end

#---------------------------------------------------------------------------------------#  
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function writeresults_pbcg_static(filename, cg_iter, rmp_obj, rmp_time, pp_time, pp_time_par, x, y, z, w, ordtime, minreducedcost)

	bigOrderArcSet = []
	for i in orders
       bigOrderArcSet = union(bigOrderArcSet, orderArcSet[i])
    end

	#Calculate necessary quantities
	totalorders = length(orders)

	obj_delay = sum((getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) 
	obj_dist = sum(sum(sum(c[a]*delta[i,a,p]*getvalue(x[i,p]) for a in orderArcSet[i]) for p in paths[i]) for i in orders) + sum(c[a]*getvalue(y[a]) for a in A_hasdriver) + sum(u[a]*getvalue(w[a]) for a in A_space) 
		
	infeasible_order_count = sum(getvalue(x[i,1]) for i in orders)
	if intersect(bigOrderArcSet, [a for a in numarcs+1:extendednumarcs]) == []
		incomplete_order_count = 0
	else
		incomplete_order_count = sum(sum(sum(delta[i,a,p] * getvalue(x[i,p]) for p in paths[i]) for a in intersect(orderArcSet[i], [a2 for a2 in numarcs+1:extendednumarcs])) for i in orders if intersect(orderArcSet[i], [a2 for a2 in numarcs+1:extendednumarcs]) != [])
	end

	delay_per_order = obj_delay/length(orders)
	max_order_delay = maximum([(getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders]) 
	totalnightsathome = sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)
	totalnightsaway = sum(length(T_off_0[d]) for d in drivers) - sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)

	totalordertrips, totalemptytrips, totaltaxitrips = 0, 0, 0
	totalordermiles, totalemptymiles, totaltaximiles = 0, 0, 0
	for i in orders, p in paths[i]
		if getvalue(x[i,p]) > 0.001
			for a in intersect(A_space, orderArcSet[i])
				if delta[i,a,p] > 0.001
					totalordertrips += delta[i,a,p] * getvalue(x[i,p])
					totalordermiles += delta[i,a,p] * c[a] *getvalue(x[i,p])
				end
			end
		end
	end
	for a in intersect(A_space, A_hasdriver)
		if getvalue(y[a]) > 0.001
			totalemptytrips += getvalue(y[a])
			totalemptymiles += c[a]*getvalue(y[a])
		end
	end

	for a in A_space
		if getvalue(w[a]) > 0.001
			totaltaxitrips += getvalue(w[a])
			totaltaximiles += u[a]*getvalue(w[a])
		end
	end

	#----------------------------WRITE ORDER TIMES TO CSV----------------------------#

	df = DataFrame(ID = [runid],
			example = [ex],
			timehorizon = [horizon], 
			timestep = [tstep], 
			numlocs = [numlocs], 
			numdrivers = [length(drivers)], 
			numtrucks = [numtrucks], 
			numarcs = [numarcs], 
			lambdaobj = [lambda],
			numorders = [totalorders], 
			pbcgiteration = [cg_iter],
			heuristicdp = [solvedpheuristically_flag_now],
			totalpbcgiteration = [""],
			minreducedcost = [minreducedcost],
			optval = [rmp_obj], 
			optval_delay = [obj_delay], 
			optval_dist = [obj_dist], 
			infeasible_order_count = [infeasible_order_count],
			incomplete_order_count = [incomplete_order_count],
			ordertrips = [totalordertrips], 
			emptytrips = [totalemptytrips], 
			taxitrips = [totaltaxitrips],
			ordermiles = [totalordermiles], 
			emptymiles = [totalemptymiles], 
			taximiles = [totaltaximiles], 
			delay_per_order = [delay_per_order], 
			max_order_delay = [max_order_delay],
			nightsathome = [totalnightsathome],
			nightsaway = [totalnightsaway],
			ip_time = [""],
			rmp_time = [rmp_time],
			pp_time = [pp_time],
			pptime_parallel = [pp_time_par],
			path_count = [sum(length(paths[i]) for i in orders)]
           )

	if cg_iter == 1
		CSV.write(filename, df)
	else
		CSV.write(filename, df, append=true)
	end

end

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function writeresults_ip_pbcg_static(filename, comptimesfilename, total_cg_iter, ip_obj, ip_time, cg_rmptimes, cg_pptimes, cg_pptimes_par, x, y, z, w, ordtime)

	bigOrderArcSet = []
	for i in orders
       bigOrderArcSet = union(bigOrderArcSet, orderArcSet[i])
    end

	#Calculate necessary quantities
	totalorders = length(orders)

	obj_delay = sum((getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) 
	obj_dist = sum(sum(sum(c[a]*delta[i,a,p]*getvalue(x[i,p]) for a in orderArcSet[i]) for p in paths[i]) for i in orders) + sum(c[a]*getvalue(y[a]) for a in A_hasdriver) + sum(u[a]*getvalue(w[a]) for a in A_space) 
		
	infeasible_order_count = sum(getvalue(x[i,1]) for i in orders)
	if intersect(bigOrderArcSet, [a for a in numarcs+1:extendednumarcs]) == []
		incomplete_order_count = 0
	else
		incomplete_order_count = sum(sum(sum(delta[i,a,p] * getvalue(x[i,p]) for p in paths[i]) for a in intersect(orderArcSet[i], [a2 for a2 in numarcs+1:extendednumarcs])) for i in orders if intersect(orderArcSet[i], [a2 for a2 in numarcs+1:extendednumarcs]) != [])
	end

	delay_per_order = obj_delay/length(orders)
	max_order_delay = maximum([(getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders]) 
	totalnightsathome = sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)
	totalnightsaway = sum(length(T_off_0[d]) for d in drivers) - sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)

	totalordertrips, totalemptytrips, totaltaxitrips = 0, 0, 0
	totalordermiles, totalemptymiles, totaltaximiles = 0, 0, 0
	for i in orders, p in paths[i]
		if getvalue(x[i,p]) > 0.001
			for a in intersect(A_space, orderArcSet[i])
				if delta[i,a,p] > 0.001
					totalordertrips += delta[i,a,p]* getvalue(x[i,p])
					totalordermiles += delta[i,a,p] * c[a] *getvalue(x[i,p])
				end
			end
		end
	end
	for a in intersect(A_space, A_hasdriver)
		if getvalue(y[a]) > 0.001
			totalemptytrips += getvalue(y[a])
			totalemptymiles += c[a]*getvalue(y[a])
		end
	end

	for a in A_space
		if getvalue(w[a]) > 0.001
			totaltaxitrips += getvalue(w[a])
			totaltaximiles += u[a]*getvalue(w[a])
		end
	end

	rmp_time = sum(cg_rmptimes)
	cg_pp_time = sum(cg_pptimes)
	cg_pptime_parallel = sum(cg_pptimes_par)

	#----------------------------WRITE ORDER TIMES TO CSV----------------------------#

	df = DataFrame(ID = [runid],
			example = [ex],
			timehorizon = [horizon], 
			timestep = [tstep], 
			numlocs = [numlocs], 
			numdrivers = [length(drivers)], 
			numtrucks = [numtrucks], 
			numarcs = [numarcs], 
			lambdaobj = [lambda],
			numorders = [length(orders)], 
			abcgiteration = ["IP"],
			heuristicdp = [solvedpheuristically_flag_now],
			totalabcgiteration = [total_cg_iter],
			minreducedcost = [""],
			optval = [ip_obj], 
			optval_delay = [obj_delay], 
			optval_dist = [obj_dist], 
			infeasible_order_count = [infeasible_order_count],
			incomplete_order_count = [incomplete_order_count],
			ordertrips = [totalordertrips], 
			emptytrips = [totalemptytrips], 
			taxitrips = [totaltaxitrips],
			ordermiles = [totalordermiles], 
			emptymiles = [totalemptymiles], 
			taximiles = [totaltaximiles], 
			delay_per_order = [delay_per_order], 
			max_order_delay = [max_order_delay],
			nightsathome = [totalnightsathome],
			nightsaway = [totalnightsaway],
			ip_time = [ip_time],
			cg_rmp_time = [rmp_time],
			cg_pp_time = [cg_pp_time],
			cg_pptime_parallel = [cg_pptime_parallel],
			path_count = [sum(length(paths[i]) for i in orders)]
           )

	CSV.write(filename, df, append=true)

	#--------------------------------------------------------------------------------#

	df2 = DataFrame(ID = [runid],
			example = [ex],
			numorders = [length(orders)], 
			numlocs = [numlocs], 
			numdrivers = [length(drivers)], 
			numtrucks = [numtrucks], 
			timehorizon = [horizon], 
			timestep = [tstep], 
			lambdaobj = [lambda],
			objval = [ip_obj],
			delayobj = [obj_delay],
			milesobj = [obj_dist],
			totaltime = [ip_time + rmp_time + cg_pptime_parallel],
			path_count = [sum(length(paths[i]) for i in orders)]
           )

	CSV.write(comptimesfilename, df2)

end

#---------------------------------------------------------------------------------------#  
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function writeresults_fullip(filename, comptimesfilename, ip_obj, totaltime_ip, orderArcSet, x, y, z, w, ordtime, currtime)

	bigOrderArcSet = []
	for i in orders
       bigOrderArcSet = union(bigOrderArcSet, orderArcSet[i])
    end

	#Calculate necessary quantities
	totalorders = length(orders)

	obj_delay = sum((getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders) 
	obj_dist = sum(sum(c[a]*getvalue(x[i,a]) for a in orderArcSet[i]) for i in orders) + sum(c[a]*getvalue(y[a]) for a in A_hasdriver) + sum(u[a]*getvalue(w[a]) for a in A_space) 
	infeasible_order_count = sum(getvalue(x[i,dummyarc]) for i in orders)
	if intersect(bigOrderArcSet, [a for a in numarcs+1:extendednumarcs]) == []
		incomplete_order_count = 0
	else
		incomplete_order_count = sum(sum(getvalue(x[i,a]) for a in intersect(orderArcSet[i], [a for a in numarcs+1:extendednumarcs])) for i in orders if intersect(orderArcSet[i], [a for a in numarcs+1:extendednumarcs]) != [])
	end

	delay_per_order = obj_delay/length(orders)
	max_order_delay = maximum([(getvalue(ordtime[i]) - shortesttriptimes[i])/shortesttriptimes[i] for i in orders]) 

	totalnightsathome = sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)
	totalnightsaway = sum(length(T_off_0[d]) for d in drivers) - sum(sum(getvalue(z[d,arcs[nodes[(driverHomeLocs[d],t2)], nodes[(driverHomeLocs[d],t2+tstep)]]]) for t2 in T_off_0[d]) for d in drivers)

	totalordertrips, totalemptytrips, totaltaxitrips = 0, 0, 0
	totalordermiles, totalemptymiles, totaltaximiles = 0, 0, 0
	for i in orders, a in intersect(A_space, orderArcSet[i])
		if getvalue(x[i,a]) > 0.001
			totalordertrips += 1
			totalordermiles += c[a]
		end
	end
	for a in intersect(A_space, A_hasdriver)
		if getvalue(y[a]) > 0.001
			totalemptytrips += getvalue(y[a])
			totalemptymiles += c[a]*getvalue(y[a])
		end
	end

	for a in A_space
		if getvalue(w[a]) > 0.001
			totaltaxitrips += getvalue(w[a])
			totaltaximiles += u[a]*getvalue(w[a])
		end
	end

	#----------------------------WRITE ORDER TIMES TO CSV----------------------------#

	df = DataFrame(ID = [runid],
			example = [ex],
			timehorizon = [horizon], 
			timestep = [tstep], 
			numlocs = [numlocs], 
			numdrivers = [length(drivers)], 
			numtrucks = [numtrucks], 
			numarcs = [numarcs], 
			lambdaobj = [lambda],
			numorders = [length(orders)], 
			iteration = ["IP"],
			optval = [ip_obj], 
			optval_delay = [obj_delay], 
			optval_dist = [obj_dist], 
			infeasible_order_count = [infeasible_order_count],
			incomplete_order_count = [incomplete_order_count],
			ordertrips = [totalordertrips], 
			emptytrips = [totalemptytrips], 
			taxitrips = [totaltaxitrips],
			ordermiles = [totalordermiles], 
			emptymiles = [totalemptymiles], 
			taximiles = [totaltaximiles], 
			delay_per_order = [delay_per_order], 
			max_order_delay = [max_order_delay],
			nightsathome = [totalnightsathome],
			nightsaway = [totalnightsaway],
			ip_time = [totaltime_ip],
			order_arc_count = [sum(length(orderArcSet[i]) for i in orders)]
           )

	CSV.write(filename, df) #, append=true)

	#--------------------------------------------------------------------------------#

	df2 = DataFrame(ID = [runid],
			example = [ex],
			numorders = [length(orders)], 
			numlocs = [numlocs], 
			numdrivers = [length(drivers)], 
			numtrucks = [numtrucks], 
			timehorizon = [horizon], 
			timestep = [tstep], 
			lambdaobj = [lambda],
			objval = [ip_obj],
			delayobj = [obj_delay],
			milesobj = [obj_dist],
			totaltime = [totaltime_ip],
			order_arc_count = [sum(length(orderArcSet[i]) for i in orders)]
           )

	CSV.write(comptimesfilename, df2)

end

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function savefullsolution_static(fullxsolutionfilename, fullysolutionfilename, fullzsolutionfilename, x, y, z, orderArcSet)

	xlist, zlist = [], []
	for i in orders, a in setdiff(orderArcSet[i], dummyarc)
		push!(xlist, (i,a))
	end
	for d in drivers, a in homeArcSet[d]
		push!(zlist, (d,a))
	end

	#Save x solution
	df_x = DataFrame(varname = ["x" for item in xlist], 
				   i = [item[1] for item in xlist],
				   a = [item[2] for item in xlist],
				   arcstartloc = [nodesLookup[arcLookup[item[2]][1]][1] for item in xlist],
				   arcendloc = [nodesLookup[arcLookup[item[2]][2]][1] for item in xlist],
				   arcstarttime = [nodesLookup[arcLookup[item[2]][1]][2] for item in xlist],
				   arcendtime =  [nodesLookup[arcLookup[item[2]][2]][2] for item in xlist],
				   value = [getvalue(x[item]) for item in xlist]
					)

	CSV.write(fullxsolutionfilename, df_x)

	#Save y solution
	df_y = DataFrame(varname = ["y" for item in A_hasdriver], 
				   a = [item for item in A_hasdriver],
				   arcstartloc = [nodesLookup[arcLookup[item][1]][1] for item in A_hasdriver],
				   arcendloc = [nodesLookup[arcLookup[item][2]][1] for item in A_hasdriver],
				   arcstarttime = [nodesLookup[arcLookup[item][1]][2] for item in A_hasdriver],
				   arcendtime =  [nodesLookup[arcLookup[item][2]][2] for item in A_hasdriver],
				   value = [getvalue(y[item]) for item in A_hasdriver]
					)

	CSV.write(fullysolutionfilename, df_y)

	#Save z solution
	df_z = DataFrame(varname = ["z" for item in zlist], 
				   i = [item[1] for item in zlist],
				   a = [item[2] for item in zlist],
				   arcstartloc = [nodesLookup[arcLookup[item[2]][1]][1] for item in zlist],
				   arcendloc = [nodesLookup[arcLookup[item[2]][2]][1] for item in zlist],
				   arcstarttime = [nodesLookup[arcLookup[item[2]][1]][2] for item in zlist],
				   arcendtime =  [nodesLookup[arcLookup[item[2]][2]][2] for item in zlist],
				   value = [getvalue(z[item]) for item in zlist]
					)

	CSV.write(fullzsolutionfilename, df_z)

end

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function savefullsolution_pbcg(fullxsolutionfilename, fullysolutionfilename, fullzsolutionfilename, x, y, z)

	xlist, zlist = [], []
	for i in orders, p in paths[i], a in setdiff(orderArcSet[i], dummyarc)
		if delta[i,a,p] > 0.01
			push!(xlist, (i,a,p))
		end
	end
	for d in drivers, a in homeArcSet[d]
		push!(zlist, (d,a))
	end

	#Save x solution
	df_x = DataFrame(varname = ["x" for item in xlist], 
				   i = [item[1] for item in xlist],
				   a = [item[2] for item in xlist],
				   arcstartloc = [nodesLookup[arcLookup[item[2]][1]][1] for item in xlist],
				   arcendloc = [nodesLookup[arcLookup[item[2]][2]][1] for item in xlist],
				   arcstarttime = [nodesLookup[arcLookup[item[2]][1]][2] for item in xlist],
				   arcendtime =  [nodesLookup[arcLookup[item[2]][2]][2] for item in xlist],
				   value = [getvalue(x[item[1], item[3]]) for item in xlist]
					)

	CSV.write(fullxsolutionfilename, df_x)

	#Save y solution
	df_y = DataFrame(varname = ["y" for item in A_hasdriver], 
				   a = [item for item in A_hasdriver],
				   arcstartloc = [nodesLookup[arcLookup[item][1]][1] for item in A_hasdriver],
				   arcendloc = [nodesLookup[arcLookup[item][2]][1] for item in A_hasdriver],
				   arcstarttime = [nodesLookup[arcLookup[item][1]][2] for item in A_hasdriver],
				   arcendtime =  [nodesLookup[arcLookup[item][2]][2] for item in A_hasdriver],
				   value = [getvalue(y[item]) for item in A_hasdriver]
					)

	CSV.write(fullysolutionfilename, df_y)

	#Save z solution
	df_z = DataFrame(varname = ["z" for item in zlist], 
				   i = [item[1] for item in zlist],
				   a = [item[2] for item in zlist],
				   arcstartloc = [nodesLookup[arcLookup[item[2]][1]][1] for item in zlist],
				   arcendloc = [nodesLookup[arcLookup[item[2]][2]][1] for item in zlist],
				   arcstarttime = [nodesLookup[arcLookup[item[2]][1]][2] for item in zlist],
				   arcendtime =  [nodesLookup[arcLookup[item[2]][2]][2] for item in zlist],
				   value = [getvalue(z[item]) for item in zlist]
					)

	CSV.write(fullzsolutionfilename, df_z)

end

#---------------------------------------------------------------------------------------#  
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function updatepastsegments_static(x, y, z, w)

	pastordersegments, pastdriversegments, pastdriversegments_space, pastemptysegments, pasttaxisegments = [], [], [], [], []

	#Find the set of arcs that are locked in (i.e. do not include dummy or extended arcs)
	lockedarcs = [a for a in 1:numarcs]

	#====================================================#

	#Add order segments
	for i in orders, a in intersect(lockedarcs, setdiff(orderArcSet[i], dummyarc))
		if getvalue(x[i,a]) > 0.001
			#Add new current segment = (i, starttime, endtime, startloc, endloc)
			arcstartloc, arcendloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
			push!(pastordersegments, (i, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], arcstartloc, arcendloc))
		end
	end

	#Add driver segments
	for d in drivers, a in intersect(lockedarcs, homeArcSet[d])
		if getvalue(z[d,a]) > 0.001
			#Add new current segment = (dd, starttime, endtime, startloc, endloc)
			arcstartloc, arcendloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
			push!(pastdriversegments, (d, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], arcstartloc, arcendloc))
			if a in A_space
				push!(pastdriversegments_space, (d, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
			end
		end
	end

	#Add empty segments
	for a in intersect(lockedarcs, A_hasdriver_space)
		if getvalue(y[a]) > 0.001
			#Add new current segment = (numtrucks, starttime, endtime, startloc, endloc)
			push!(pastemptysegments, (getvalue(y[a]), nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
		end
	end

	#Add taxi segments
	for a in intersect(lockedarcs, A_space)
		if getvalue(w[a]) > 0.001
			#Add new current segment = (numtrucks, starttime, endtime, startloc, endloc)
			push!(pasttaxisegments, (getvalue(w[a]), nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
		end
	end

	return pastordersegments, pastdriversegments, pastdriversegments_space, pastemptysegments, pasttaxisegments

end

#---------------------------------------------------------------------------------------#  
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function updatepastextendedsegments_static(x)

	pastorderextendedsegments = []

	#Find the set of extended arcs
	lockedarcs = [a for a in numarcs+1:extendednumarcs]

	#====================================================#

	#Add order segments
	for i in orders, a in intersect(lockedarcs, setdiff(orderArcSet[i], dummyarc))
		if getvalue(x[i,a]) > 0.001
			#Add new current segment = (i, starttime, endtime, startloc, endloc)
			arcstartloc, arcendloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
			push!(pastordersegments, (i, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], arcstartloc, arcendloc))
		end
	end

	return pastorderextendedsegments

end

#---------------------------------------------------------------------------------------#  
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function updatepastsegments_pbcg_static(x, y, z, w)

	pastordersegments, pastdriversegments, pastdriversegments_space, pastemptysegments, pasttaxisegments = [], [], [], [], []

	#Find the set of arcs that are locked in (i.e. do not include dummy or extended arcs)
	lockedarcs = [a for a in 1:numarcs]

	#====================================================#

	#Add order segments
	for i in orders, p in paths[i]
		if getvalue(x[i,p]) > 0.001
			for a in intersect(lockedarcs, setdiff(orderArcSet[i], dummyarc))
				if delta[i,a,p] > 0.001
					#Add new current segment = (i, starttime, endtime, startloc, endloc)
					arcstartloc, arcendloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
					push!(pastordersegments, (i, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], arcstartloc, arcendloc))
				end
			end
		end
	end

	#Add driver segments
	for d in drivers, a in intersect(lockedarcs, homeArcSet[d])
		if getvalue(z[d,a]) > 0.001
			#Add new current segment = (dd, starttime, endtime, startloc, endloc)
			arcstartloc, arcendloc = nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1]
			push!(pastdriversegments, (d, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], arcstartloc, arcendloc))
			if a in A_space
				push!(pastdriversegments_space, (d, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
			end
		end
	end

	#Add empty segments
	for a in intersect(lockedarcs, A_hasdriver_space)
		if getvalue(y[a]) > 0.001
			#Add new current segment = (numtrucks, starttime, endtime, startloc, endloc)
			push!(pastemptysegments, (getvalue(y[a]), nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
		end
	end

	#Add taxi segments
	for a in intersect(lockedarcs, A_space)
		if getvalue(w[a]) > 0.001
			#Add new current segment = (numtrucks, starttime, endtime, startloc, endloc)
			push!(pasttaxisegments, (getvalue(w[a]), nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ))
		end
	end

	return pastordersegments, pastdriversegments, pastdriversegments_space, pastemptysegments, pasttaxisegments

end

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#  

function writeorderoutcomesfile_static(orderoutcomesfilename, solutionmethod)

	df = DataFrame(ID = [runid for i in 1:highestorderindex],
		example = [ex for i in 1:highestorderindex],
		method = [solutionmethod for i in 1:highestorderindex],
		lambda = [lambda for i in 1:highestorderindex],
		orderid = [i for i in 1:highestorderindex],
		originloc = [originloc[i] for i in 1:highestorderindex],
		destloc = [destloc[i] for i in 1:highestorderindex],
		orderdist = ordermilesoutcomes[1:highestorderindex],
		orderpenaltydist = ordermilespenalty[1:highestorderindex], 
		orderdeliverytime = orderdelayoutcomes[1:highestorderindex],
		shortesttriptime = shortesttriptimes[1:highestorderindex]
       )

	CSV.write(orderoutcomesfilename, df)

end


#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#  
#---------------------------------------------------------------------------------------#

function writecorridorcounterfile_static(solutionmethod)

	segmentcounterDict = Dict()
	for pa in prearcs
		o, d, rddtime, rawtime, miles = pa
		segmentcounterDict[o, d] = [rddtime, rawtime, miles, 0, 0, 0, 0]
	end
	for sgmt in pastordersegments
		o, d = sgmt[4], sgmt[5]
		if o != d
			try
				segmentcounterDict[o, d][4] += 1
				segmentcounterDict[o, d][5] += 1
			catch
				1+1
			end
		end
	end
	for sgmt in pastemptysegments
		ttltrucks, o, d = sgmt[1], sgmt[4], sgmt[5]
		if o != d
			try
				segmentcounterDict[o, d][4] += ttltrucks
				segmentcounterDict[o, d][6] += ttltrucks
			catch
				1+1
			end
		end
	end
	for sgmt in pasttaxisegments
		ttltaxis, o, d = sgmt[1], sgmt[4], sgmt[5]
		if o != d
			try
				segmentcounterDict[o, d][4] += ttltaxis
				segmentcounterDict[o, d][7] += ttltaxis
			catch
				1+1
			end
		end
	end

	df = DataFrame(ID = [runid for a in 1:length(segmentcounterDict)],
				example = [ex for a in 1:length(segmentcounterDict)],
				lambda = [lambda for a in 1:length(segmentcounterDict)],
				method = [solutionmethod for a in 1:length(segmentcounterDict)],
				originloc = [item[1][1] for item in segmentcounterDict],
				destloc = [item[1][2] for item in segmentcounterDict],
				arctime_rdd = [item[2][1] for item in segmentcounterDict],
				arctime_raw = [item[2][2] for item in segmentcounterDict],
				arcdist = [item[2][3] for item in segmentcounterDict],
				totaltrips = [item[2][4] for item in segmentcounterDict],
				ordertrips = [item[2][5] for item in segmentcounterDict],
				emptytrips = [item[2][6] for item in segmentcounterDict],
				taxitrips = [item[2][7] for item in segmentcounterDict]
	           )

	CSV.write(corridorfilename, df)

end

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function staticviz(x, y, z, w)

	orderarcsegments = []

	#Add order arc segments (arcs generated through column generation)
	#for i in orders, a in setdiff(orderArcSet[i], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
	for i in orders, a in setdiff([a2 for a2 in 1:numarcs], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
		#Add new current segment = (i, starttime, endtime, startloc, endloc)
		push!(orderarcsegments, (i, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ) )
	end

	homearcsegments = []
	#Add driver arc segments
	for d in drivers, a in setdiff(homeArcSet[d], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
		#Add new current segment = (d, starttime, endtime, startloc, endloc)
		push!(homearcsegments, (d, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ) )
	end

	#====================================================#

	#Visualize
	if maketimespacevizfiles == 1
		timespaceviz_online(string(vizfoldername, "/TimeSpaceNetworkMaps/", vizfilename, ".png"), showdrivers_flag, 0, nodesLookup, nodes, timeperiods, [], [], [], [], [])
	end

	if makeadvancedvizfiles == 1
		for i in 1:highestorderindex
			timespaceviz_online_byorder(string(vizfoldername, "/OrderMaps/", vizfilename, "order", i, ".png"), i, 0, nodesLookup, nodes, timeperiods, [], orderarcsegments)
		end

		for d in drivers
			timespaceviz_online_bydriver(string(vizfoldername, "/DriverMaps/", vizfilename, "driver", d, ".png"), d, 0, nodesLookup, nodes, timeperiods, [], homearcsegments)
		end
	end

end

#---------------------------------------------------------------------------------------#

function staticviz_abcgiter(x, y, z, w, cg_iter, arcredcosts)

	currentordersegments = []

	#Add order arc segments (arcs generated through column generation)
	for i in orders, a in setdiff(orderArcSet[i], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
	#for i in orders, a in setdiff([a2 for a2 in 1:numarcs], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
		#Add new current segment = (i, starttime, endtime, startloc, endloc)
		if getvalue(x[i,a]) > 0.001
			push!(currentordersegments, (i, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ) )
		end
	end

	#====================================================#

	orderarcsegments = []

	#Add order arc segments (arcs generated through column generation)
	for i in orders, a in setdiff(orderArcSet[i], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
	#for i in orders, a in setdiff([a2 for a2 in 1:numarcs], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
		#Add new current segment = (i, starttime, endtime, startloc, endloc)
		push!(orderarcsegments, (i, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ) )
	end

	homearcsegments = []
	#Add driver arc segments
	for d in drivers, a in setdiff(homeArcSet[d], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
		#Add new current segment = (d, starttime, endtime, startloc, endloc)
		push!(homearcsegments, (d, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ) )
	end

	#====================================================#

	negrecarcs = []

	for i in orders, a in setdiff(orderArcSet_viz[i], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
		if arcredcosts[i,a] < -0.001
			push!(negrecarcs, (i, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ) )
		end
	end

	#====================================================#

	#Visualize
	if makeadvancedvizfiles == 1
		for i in 1:highestorderindex
			timespaceviz_online_byorder(string(vizfoldername, "/OrderMaps/", vizfilename, "order", i, "_iter", cg_iter, ".png"), i, 0, nodesLookup, nodes, timeperiods, currentordersegments, orderarcsegments)
		end

		for i in 1:highestorderindex
			timespaceviz_negrc_byorder(string(vizfoldername, "/OrderMaps/", vizfilename, "negrc_order", i, "_iter", cg_iter, ".png"), i, 0, nodesLookup, nodes, timeperiods, currentordersegments, orderarcsegments, negrecarcs)
		end
	end

end

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#

function staticviz_pbcg(x, y, z, w)

	orderarcsegments = []

	#Add order arc segments (arcs generated through column generation)
	for i in orders, p in paths[i], a in setdiff(orderArcSet[i], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
		#Add new current segment = (i, starttime, endtime, startloc, endloc)
		if delta[i,a,p] > 0.001
			push!(orderarcsegments, (i, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ) )
		end
	end

	homearcsegments = []
	#Add driver arc segments
	for d in drivers, a in setdiff(homeArcSet[d], union(dummyarc, [a for a in numarcs+1:extendednumarcs]))
		#Add new current segment = (d, starttime, endtime, startloc, endloc)
		push!(homearcsegments, (d, nodesLookup[arcLookup[a][1]][2], nodesLookup[arcLookup[a][2]][2], nodesLookup[arcLookup[a][1]][1], nodesLookup[arcLookup[a][2]][1] ) )
	end

	#====================================================#

	#Visualize
	if maketimespacevizfiles == 1
		timespaceviz_online(string(vizfoldername, "/Time Space Network Maps/", vizfilename, ".png"), showdrivers_flag, 0, nodesLookup, nodes, timeperiods, [], [], [], [], [])
	end

	if makeadvancedvizfiles == 1
		for i in 1:highestorderindex
			timespaceviz_online_byorder(string(vizfoldername, "/Order Maps/", vizfilename, "order", i, ".png"), i, 0, nodesLookup, nodes, timeperiods, [], orderarcsegments)
		end

		for d in drivers
			timespaceviz_online_bydriver(string(vizfoldername, "/Driver Maps/", vizfilename, "driver", d, ".png"), d, 0, nodesLookup, nodes, timeperiods, [], homearcsegments)
		end
	end

end

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#