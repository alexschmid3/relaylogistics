
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
			method = [solutionmethod],
			variablefixingthreshold = [variablefixingthreshold],
			iteration = [iteration],
			objective = [obj],
			smptime = [timeslist.mp],
			pptime = [timeslist.pp],
			pptime_par = [timeslist.pppar],
			iptime = [timeslist.ip],
            totalarcs = [totalarcs]
		)

	if appendflag == 1
		CSV.write(resultsfilename, df, append=true)
	else
		CSV.write(resultsfilename, df)
	end

end
