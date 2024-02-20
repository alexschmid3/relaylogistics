
function calcdrivermetrics(z)

	util = sum(sum(fragworkinghours[driverHomeLocs[d],drivershift[d],f] * value(z[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers) / sum(maxweeklydriverhours for d in drivers)
	nightsaway = sum(sum(fragmentnightsaway[driverHomeLocs[d],drivershift[d],f] * value(z[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers) 
	
	hoursworkedperdriver = [sum(fragworkinghours[driverHomeLocs[d],drivershift[d],f] * value(z[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers] 
	driversunused = length([d for d in hoursworkedperdriver if d==0])

	return util, nightsaway, driversunused

end

#----------------------------------------------------------------------------------------#

function calcordermetrics(x)

	return 0

end

#----------------------------------------------------------------------------------------#

function writeresultsforearlytests(resultsfilename, appendflag, iteration, obj, timeslist, totalarcs, x, z)

	if formulation == "heterogeneous"
		util, nightsaway, driversunused = calcdrivermetrics(z)
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
			driversunused = [driversunused]
		)

	if appendflag == 1
		CSV.write(resultsfilename, df, append=true)
	else
		CSV.write(resultsfilename, df)
	end

end

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
