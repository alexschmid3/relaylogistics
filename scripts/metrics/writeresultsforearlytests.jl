
function writeresultsforearlytests(resultsfilename, appendflag, iteration, obj, timeslist, totalarcs)

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
			cuttime = [timeslist.cut]
		)

	if appendflag == 1
		CSV.write(resultsfilename, df, append=true)
	else
		CSV.write(resultsfilename, df)
	end

end
