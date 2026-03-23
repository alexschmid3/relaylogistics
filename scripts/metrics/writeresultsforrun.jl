
function calcdrivermetrics(z)

	hoursworked = [sum(fragworkinghours[driverHomeLocs[d],drivershift[d],f] * value(z[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers] 
	util = sum(hoursworked[d] for d in drivers) / sum(maxweeklydriverhours for d in drivers)
	nightsaway = sum(sum(fragmentnightsaway[driverHomeLocs[d],drivershift[d],f] * value(z[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers) 
	
	hoursworkedperdriver = [sum(fragworkinghours[driverHomeLocs[d],drivershift[d],f] * value(z[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers] 
	driversunused = length([d for d in hoursworkedperdriver if d==0])

	if driverhourshistogram_flag == 1
		df = DataFrame(experiment_id = [experiment_id for d in drivers], formulation = [formulation for d in drivers], driver = [d for d in drivers], shift = drivershift, homeloc = [driverHomeLocs[d] for d in drivers], hours = hoursworked)
		CSV.write(driverhourshistogramfilename, df)
	end

	return util, nightsaway, driversunused

end

#----------------------------------------------------------------------------------------#

function calcdrivermetrics_homogeneous(z)

	hoursworkedperdriver, listofdrivers, listofhomelocs, listofshifts = [], [], [], []
	for l in 1:numlocs, s in 1:numshifts
		#Solve little MIP to equitably distribute drivers among fragments
		possiblefragments = [f for f in 1:numfragments[l,s] if value(z[l,s,f]) > 1e-4]

		ip = Model(Gurobi.Optimizer)
		set_optimizer_attribute(ip, "TimeLimit", 30)
		set_optimizer_attribute(ip, "OutputFlag", 0)
		set_optimizer_attribute(ip, "MIPGap", 0.01)

		#Variables
		@variable(ip, z_new[d = driversets[l,s], f = 1:numfragments[l,s]] >= 0, Int)
		@variable(ip, maxhours >= 0)

		#Objective
		@objective(ip, Min, maxhours)

		#Driver constraints
		@constraint(ip, driverStartingLocs[d in driversets[l,s]], sum(sum(z_new[d,f] for f in F_plus_ls[l,s,n]) for n in driverSetStartNodes[l,s]) == 1)
		@constraint(ip, driverFlowBalance[d in driversets[l,s], n in N_flow_ls[l,s]], sum(z_new[d,f] for f in F_minus_ls[l,s,n]) - sum(z_new[d,f] for f in F_plus_ls[l,s,n]) == 0)
		@constraint(ip, driverMaxHours[d in driversets[l,s]], sum(fragworkinghours[l,s,f] * z_new[d,f] for f in possiblefragments) <= maxhours)
		for f in setdiff(1:numfragments[l,s], possiblefragments), d in driversets[l,s]
			@constraint(ip, z_new[d,f] == 0)
		end
		#@constraint(ip, maxhours <= maxweeklydriverhours)
		@constraint(ip, coverallfragments[f in 1:numfragments[l,s]], sum(z_new[d,f] for d in driversets[l,s]) == round(value(z[l,s,f]), digits = 0))

		#Solve to assign drivers to fragments
		optimize!(ip)

		#Record orders for each driver
		for d in driversets[l,s]
			push!(hoursworkedperdriver, sum(fragworkinghours[l,s,f] * value(z_new[d,f]) for f in possiblefragments) )
			push!(listofdrivers, d)
			push!(listofhomelocs, l)
			push!(listofshifts, s)
		end
	end
	driversunused = length([d for d in hoursworkedperdriver if d==0])

	util = sum(hoursworkedperdriver) / sum(shiftlength*horizon/24 for d in drivers)
	nightsaway = sum(sum(sum(fragmentnightsaway[l,s,f] * value(z[l,s,f]) for f in 1:numfragments[l,s]) for l in 1:numlocs if length(driversets[l,s]) > 1e-4) for s in 1:numshifts)

	if driverhourshistogram_flag == 1
		df = DataFrame(experiment_id = [experiment_id for d in drivers], formulation = [formulation for d in drivers], driver = listofdrivers, shift = listofshifts, homeloc = listofhomelocs, hours = hoursworkedperdriver)
		CSV.write(driverhourshistogramfilename, df)
	end

	return util, nightsaway, driversunused

end

#----------------------------------------------------------------------------------------#

function calcordermetrics(x)

	return 0

end

#----------------------------------------------------------------------------------------#

function writeresultsforrun(resultsfilename, appendflag, iteration, obj, timeslist, totalarcs, x, z)

	if iteration == "IP"
		if (formulation == "heterogeneous") & !(solutionmethod == "arcip")
			util, nightsaway, driversunused = calcdrivermetrics(z)
		else
			util, nightsaway, driversunused = calcdrivermetrics_homogeneous(z)
		end
	else
		util, nightsaway, driversunused = 0, 0, 0
	end
	
	df = DataFrame(experiment_id = [experiment_id],
			instance = [ex],
			lambda_delay = [lambda],
			lambda_drvrhrs = [lambda2],
			horizon = [horizon],
			tstep = [tstep],
			week = [weekstart],
			numlocs = [numlocs],
			numorders = [length(orders)],
			numdrivers = [length(drivers)],
			maxweeklydriverhours = [maxweeklydriverhours],
			method = [solutionmethod],
			variablefixingthreshold = [variablefixingthreshold],
			varsettingiterations = [varsettingiterations],
			strongreducedcosts = [strengthenedreducedcost_flag],
			columnmemory = [columnmemorylength],
			deletioniterationpercent = [postmagcolumndeletioniterationpercent],
			deletionthreshold = [postmagcolumndeletionthreshold],
			cuttype = [knapsackcuttype],
			iteration = [iteration],
			objective = [obj],
			smptime = [timeslist.mp],
			pptime = [timeslist.pp],
			pptime_par = [timeslist.pppar],
			iptime = [timeslist.ip],
            totalarcs = [totalarcs],
			cuttime = [timeslist.cut],
			driverutil = [util],
			drivernightsaway = [nightsaway],
			driversunused = [driversunused],
			fulltime = [timeslist.full]
		)

	if appendflag == 1
		CSV.write(resultsfilename, df, append=true)
	else
		CSV.write(resultsfilename, df)
	end

end

#----------------------------------------------------------------------------------------#

function writeresultsforrun_deadlines(resultsfilename, appendflag, iteration, obj, timeslist, totalarcs, x, z, totalmiles, totaldelay, totalordertime, totalemptymiles, totalrepomiles, totalshortestpathmiles, totalordermiles, totalpenaltymiles, timelaboring, timeasinventory)

	if iteration == "IP"
		if (formulation == "heterogeneous") & !(solutionmethod == "arcip")
			util, nightsaway, driversunused = calcdrivermetrics(z)
		else
			util, nightsaway, driversunused = calcdrivermetrics_homogeneous(z)
		end
	else
		util, nightsaway, driversunused = 0, 0, 0
	end
	
	df = DataFrame(experiment_id = [experiment_id],
			instance = [ex],
			lambda_delay = [lambda],
			lambda_drvrhrs = [lambda2],
			horizon = [horizon],
			tstep = [tstep],
			week = [weekstart],
			numlocs = [numlocs],
			numorders = [length(orders)],
			numdrivers = [length(drivers)],
			maxweeklydriverhours = [maxweeklydriverhours],
			method = [solutionmethod],
			variablefixingthreshold = [variablefixingthreshold],
			varsettingiterations = [varsettingiterations],
			strongreducedcosts = [strengthenedreducedcost_flag],
			columnmemory = [columnmemorylength],
			deletioniterationpercent = [postmagcolumndeletioniterationpercent],
			deletionthreshold = [postmagcolumndeletionthreshold],
			cuttype = [knapsackcuttype],
			iteration = [iteration],
			objective = [obj],
			smptime = [timeslist.mp],
			pptime = [timeslist.pp],
			pptime_par = [timeslist.pppar],
			iptime = [timeslist.ip],
            totalarcs = [totalarcs],
			cuttime = [timeslist.cut],
			driverutil = [util],
			drivernightsaway = [nightsaway],
			driversunused = [driversunused],
			fulltime = [timeslist.full],
			totalmiles = [totalmiles], 
			totaldelay = [totaldelay], 
			totalordertime = [totalordertime], 
			totalemptymiles = [totalemptymiles], 
			totalrepomiles = [totalrepomiles],
			deadline_sp = [deadlineasmultipleofshortestpath],
			totalshortestpathmiles = [totalshortestpathmiles], 
			totalordermiles = [totalordermiles], 
			totalpenaltymiles = [totalpenaltymiles], 
			timelaboring = [timelaboring], 
			timeasinventory = [timeasinventory],
			percentnightshift = [percentnightshift], 
			laborcost_delta = [laborcost_delta], 
			driverinventorycost_theta = [driverinventorycost_theta] 
		)

	if appendflag == 1
		CSV.write(resultsfilename, df, append=true)
	else
		CSV.write(resultsfilename, df)
	end

end

#----------------------------------------------------------------------------------------#

function writedriverstats(filename, z)

	mylocs, myshifts, mydrivers, myunused, myhours, myutil, mynights = [], [], [], [], [], [], []
	for l in 1:numlocs, s in 1:numshifts

		push!(mylocs, l)
		push!(myshifts, s)
		push!(mydrivers, length(driversets[l,s]))
			
		if driversets[l,s] != []
			util = sum(sum(fragworkinghours[l,s,f] * value(z[d,f]) for f in 1:numfragments[l,s]) for d in driversets[l,s]) / sum(maxweeklydriverhours for d in driversets[l,s])
			nightsaway = sum(sum(fragmentnightsaway[l,s,f] * value(z[d,f]) for f in 1:numfragments[l,s]) for d in driversets[l,s]) 
			hoursworkedperdriver = [sum(fragworkinghours[l,s,f] * value(z[d,f]) for f in 1:numfragments[l,drivershift[d]]) for d in driversets[l,s]] 
			driversunused = length([d for d in hoursworkedperdriver if d==0])		
			push!(myunused, driversunused)
			push!(myhours, sum(hoursworkedperdriver))
			push!(myutil, util)
			push!(mynights, nightsaway)
		else
			push!(myunused, 0)
			push!(myhours, 0)
			push!(myutil, 1)
			push!(mynights, 0)
		end
	end

	df = DataFrame(location=mylocs,
					shift=myshifts,
					drivers=mydrivers,
					unused=myunused,
					totalhoursworked=myhours,
					utilization=myutil,
					nightsaway=mynights)
					
	CSV.write(filename, df)

end
