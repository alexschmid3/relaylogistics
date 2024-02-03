
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

	util, nightsaway, driversunused = calcdrivermetrics(z)
	
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
