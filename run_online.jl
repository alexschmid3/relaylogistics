
using JuMP, Gurobi, Random, CSV, DataFrames, Statistics, Dates, SparseArrays 

#-----------------------------------LOAD OTHER FILES------------------------------------#

include("scripts/instancegeneration/readrivigodata.jl")
include("scripts/instancegeneration/shortestpath.jl")
include("scripts/instancegeneration/constraintmatrix.jl")
include("scripts/instancegeneration/completeinstance.jl")
include("scripts/journeybasedmodel/initializejourneymodel.jl")
include("scripts/metrics/writeresultsforrun.jl")
include("scripts/onlineimplementation/initializecurrentstatearcs.jl")
include("scripts/onlineimplementation/getinitialstate.jl")
include("scripts/instancegeneration/initializearcsets.jl")
include("scripts/multiarcgeneration/initializeorderarcsets.jl")

#-------------------------------------FOUR INSTANCES------------------------------------#  

inittime = time()

loclist = [20, 40, 60, 66, 66]
driverlist = [70, 300, 1600, 3316, 1600]
trucklist = [60, 230, 1200, 2495, 1200]
seedlist = [202481, 155702, 731761, 963189, 731762]
ordercapslist = [3, 6, 18, 10000, 18]  #Cap on the number of orders introduced in each 6 hour increment in the online problem (6 hrs doesn't depend on tstep or timedelta to ensure consistent instance when these parameters are adjusted)
								       #Order caps is technically an online parameter, but we use it in the static version as well to ensure the orders included in the instance are the same in the static and online version

hubdistancesfilename = "data/hubdistances.csv"
traveltimesfilename = "data/traveltimes_outliers.csv"
hubdataisbfilename = "data/hub_data_isb_connect.csv"
vntdataisbfilename = "data/vnt_data_isb_connect_clean.csv"
lhdataisbfilename = "data/lh_data_isb_connect_clean.csv"

#-----------------------------------ONLINE PARAMETERS-----------------------------------#  	

becomesavailablehours = 24*3 	
onlineweeks = 2 #expparms[experiment_id, 4]
timedelta = 24 #expparms[experiment_id, 5]
onlinetimehorizon = 24*7*onlineweeks		  	    			# Length of the online time horizon in hours (ex. run 1 week of online iterations), should be multiple of timedelta
numiterations_online = convert(Int64, onlinetimehorizon/timedelta)		
pruneorderarcs_flag = 1								# 1 = prune order arcs that are no longer usable after each online iteration, 0 = do not
maxorderdelayhours = onlinetimehorizon * 2

#----------------------------------INSTANCE PARAMETERS----------------------------------#  	

#Read experiment parameters from file
experiment_id = 1 #ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
paramsfilename = "data/onlineruns.csv"
expparms = CSV.read(paramsfilename, DataFrame)
formulation = "homogeneous" #expparms[experiment_id, 15]  # Drivers = homogeneous, heterogeneous
ex = 1 #expparms[experiment_id, 2]		
weekstart = expparms[experiment_id, 4]
horizon = expparms[experiment_id, 5] * 24
tstep = expparms[experiment_id, 6]
driverfactor = 1.7 # expparms[experiment_id, 7]
k = expparms[experiment_id, 8]   #This is the rho value
opt_gap = expparms[experiment_id, 9]
lambda = expparms[experiment_id, 10]
maxweeklydriverhours = expparms[experiment_id, 11]
lambda2 = expparms[experiment_id, 12]
println("Experiment = ", experiment_id)

#Read algorithm control parameters from file
solutionmethod = expparms[experiment_id, 3]		
variablefixing_ub = expparms[experiment_id, 18]
variablefixing_lb = expparms[experiment_id, 17]
variablefixingthreshold = (variablefixing_lb, variablefixing_ub)
varsettingiterations = expparms[experiment_id, 13]
varsettingtype = "z" 
strengthenedreducedcost_flag = expparms[experiment_id, 14]
columnmemorylength = expparms[experiment_id, 16] #Forget unused columns after this many iterations
postmagcolumndeletioniterationpercent = expparms[experiment_id, 19] 
postmagcolumndeletionthreshold = expparms[experiment_id, 20] 
knapsackcuttype = expparms[experiment_id, 21] 
symmetrybreaking_flag = expparms[experiment_id, 22] 
if knapsackcuttype > 0
	knapsackcuts_flag = 1
else
	knapsackcuts_flag = 0
end

#Transform date
weekstart = DateTime(weekstart) + Dates.Hour(8)

#Definition of the instance 
iterationordercap = ordercapslist[ex]					 
maxlocs = loclist[ex]
maxdrivers = round(driverlist[ex] / driverfactor, digits = 0)							 
numtrucks = round(trucklist[ex] / driverfactor, digits = 0)
randomseedval = seedlist[ex]
Random.seed!(randomseedval)

#Problem/optimization parameters (all hard coded values we decided on at the beginning)
shiftlength = 12									# Length of each driver shift in hours
taxicostpct = 2.0                          			# Cost to taxi along an arc as a percentage of cost to drive a truck along that arc (should be > 1.0)
roundup_flag = 1			 					 	# 0 = round down to find discretized travel times, 1 = round up
drivershifttstep = 12								# How many hours between start of driver shifts, (ex. drivershifttstep=12 means each driver's first shift starts at time 0 or time 12, drivershifttstep=6 means a driver's first shift could start at time 0, 6, 12, or 18)
tstepforordercreation = 12 							# Should be same as timedelta in equivalent online instance; used to round order available time stamps from Rivigo data (ex. round observed available time to previous 12 hour block)
inprogressdummyarc_flag = 0						# 1 = allow in progress orders to be assigned to the dummy arc, 0 = do not (Should be assigned to 0 to ensure feasibility/continuity of online iterations)
truearcfinishtime_flag = 0							# 1 = use unrounded arc travel times to assess order delivery delay (still has some bugs), 0 = use travel times rounded up to the next time step
finallegdistancepenalty = 0.40						# Distance penalty assessed for orders that finish beyond the planning horizon
finallegtimepenalty = 0.30							# Time/delay penalty assessed for orders that finish beyond the planning horizon
dummyendtime = 1000									# Dummy time assigned to the "beyond the horizon" nodes
maxnightsaway = 1
driveroffdays_flag = 0
vizflag = 0
saveconvergencedata_flag = 1

#Travel time calculation parameters
excludeoutliers_flag = 1							# 0 = include outliers in travel time calculation, 1 = exclude outliers
googlemapstraveltimes_flag = 1						# 1 = travel time between two locations is max(Avg from Rivigo data, Google Maps travel time) (<5% of arcs use google maps time), 0 = travel time is avg of Rivigo data
includesymmetricarcs_flag = 1						# 1 = if arc A-->B present in Rivigo data but not B-->A, create synthetic arc B-->A; 0 = do not include synthetic arcs (may cause feasibility issues)
traveltimefordelay_flag = 2 						# 0 = use rounded travel times for shortest path used in delay objective, 1 = use raw travel times (best for comparing across multiple tsteps), 2 = use rounded travel times, except on the final leg of the journey where raw is used 
ensureconnectivity_flag = 1

#Uniform k and ABCG + k control parameters
ktype_flag = "pct"									# "hrs" = # of hours acceptable delay, "pct" = acceptable delay as percent of shortest path time, "min24" = max(24 hrs, percent of shortest path)

#Create id for this run
runid = string("ex", ex, "_exp", experiment_id, "_", solutionmethod, "_rundate", today())

#File names					
vizfoldername = string("visualizations/static/run ", runid)
csvfoldername = string("outputs/table3/")
resultsfilename = string(csvfoldername, runid, "_output.csv")
convergencedatafilename = string(csvfoldername, "convergence_exp", runid, ".csv")

#---------------------------------IMPORT FINAL SCRIPTS----------------------------------# 

#Import algorithm functions
include("scripts/journeybasedmodel/solvejourneymodel_online.jl")
include("scripts/multiarcgeneration/multiarcgeneration_online.jl")

#Import miscellaneous functions
if vizflag == 1
	include("scripts/visualizations/timespacenetwork.jl")
end
include("scripts/directoryinitialization.jl")
include("scripts/helper.jl")

#-----------------------------------GENERATE INSTANCE-----------------------------------# 

#Initialize currentdatetime
currentdatetime = weekstart

#Create node and arc networks
timeperiods = horizon / tstep + 1
maxviztimeperiods = timeperiods 
hubCoords, hubsLookup, hubsReverseLookup, hubsTravelTimeIndex, numlocs = readlocations(hubdataisbfilename, maxlocs)
nodes, nodesLookup, N_0, N_end, numnodes = timespacentwk(numlocs, tstep, horizon)
prearcs, arcLength, arcLength_raw = readarcs(traveltimesfilename, hubdistancesfilename, tstep, numlocs, hubsTravelTimeIndex, roundup_flag, excludeoutliers_flag, hubsReverseLookup, googlemapstraveltimes_flag)
arcs, arcLookup, A_plus, A_minus, A_space, A_plus_time, A_minus_time, A_minus_space, A_plus_space, numarcs, truetraveltime = arccreation(prearcs, horizon, tstep, numnodes, nodes, numlocs)
c, u = calcobjectivecosts(hubdistancesfilename)

#Get initial state of system
currstate, includeorderidlist, drivers, driverHomeLocs, drivershift, N_flow_t, T_off_Monday8am, numshifts, originloc, destloc, orderOriginalStartLoc, orderOriginalStartTime, highestorderindex, distbetweenlocs, shortestpatharclists, traveltimebetweenlocs_rdd, traveltimebetweenlocs_raw, traveltimebetweenlocs_llr, nodesLookup, arcLookup, A_minus, A_plus, c, extendednodes, extendednumnodes, extendedarcs, extendednumarcs = getinitialstate(nodesLookup, arcLookup, A_minus, A_plus, c)

#Finish extending the time-space network
arcLookup, nodesLookup, arcfinishtime, dummyarc, allarcs = calcarcfinishtimes()

#Create initial arc sets
currarcs, currfragments, primaryarcs, extendedtimearcs, numeffshifts = initializecurrentstatearcs(currstate)

#---------------------------------------SOLVE----------------------------------------# 

for currtime in 0:timedelta:0 #timedelta*(numiterations_online-1)

	println("------------------------------------BEGIN ITERATION CURRTIME = $currtime------------------------------------")
	#=
    #Solve current instance
    if (solutionmethod == "mag") || (solutionmethod == "sag")
        mag_obj, smp, x_smp, y_smp, z_smp, w_smp, magarcs, smptime, pptime, pptime_par, totalmagarcs = multiarcgeneration!(currstate, currfragments, currarcs)    
        magip_obj, x_magip, z_magip, magip_time, magip_bound = solvejourneymodel(0, opt_gap, currstate, currarcs, currfragments, magarcs)
    end

	if (solutionmethod == "mvg") || (solutionmethod == "otdf") || (solutionmethod == "rrf")
		
		#Save most recently executed segments for visualization and reporting
		driverarcstaken = updatepastsegments_fragment(timedelta, x_ip, y_ip, z_ip, w_ip, restrictedArcSet, driversets, driversingroup)
		
		#Update whether each driver was at home or away from home during their last off hours
		updateawaylastnight_fragment(z_ip, driverarcstaken)
		
	end


    
	if (solutionmethod == "mvg") || (solutionmethod == "otdf") || (solutionmethod == "rrf")
		#Save most recently executed segments for visualization and reporting
		driverarcstaken = updatepastsegments_fragment(timedelta, x_ip, y_ip, z_ip, w_ip, restrictedArcSet, driversets, driversingroup)
		#Update whether each driver was at home or away from home during their last off hours
		updateawaylastnight_fragment(z_ip, driverarcstaken)
	elseif solutionmethod == "ptp"
		driverarcstaken = updatepastsegments_ptp(timedelta, x_ip, y_ip, z_ip, w_ip, orderArcSet_full, driversets, driversingroup)
	elseif solutionmethod != "bns"
		#Save most recently executed segments for visualization and reporting
		updatepastsegments(timedelta, x_ip, y_ip, z_ip, w_ip, restrictedArcSet)
		#Update whether each driver was at home or away from home during their last off hours
		updateawaylastnight(z_ip)
	end

	if (solutionmethod == "mvg") || (solutionmethod == "otdf") || (solutionmethod == "rrf")
		savefullsolution_online_fragment(fullxsolutionfilename, fullysolutionfilename, fullzsolutionfilename, currtime, x_ip, y_ip, z_ip, restrictedArcSet, driverarcstaken)
	elseif solutionmethod == "ptp"
		savefullsolution_online_fragment(fullxsolutionfilename, fullysolutionfilename, fullzsolutionfilename, currtime, x_ip, y_ip, z_ip, orderArcSet_full, driverarcstaken)
		#savefullsolution_online(fullxsolutionfilename, fullysolutionfilename, fullzsolutionfilename, currtime, x_ip, y_ip, z_ip, orderArcSet_full)
	elseif solutionmethod != "bns"
		savefullsolution_online(fullxsolutionfilename, fullysolutionfilename, fullzsolutionfilename, currtime, x_ip, y_ip, z_ip, restrictedArcSet)
	end

	#Update datetime
	global currentdatetime = currentdatetime + Dates.Hour(timedelta)

	#Reset index sets
	global driversintransit = []
	global trucksintransit = []
	global T_off = []
	global T_off_0 = Dict()
	global T_off_constr = Dict()
	global A_hasdriver = []
	global yupperbound = []
	global A_hasdriver_space = []

	#Iterate forward by timedelta 
	if solutionmethod == "bns"
		updatedrivers_bns(timedelta, currentdatetime, weekstart, T_off_Monday8am, driverassignments)
		updateawaylastnight_bns(driverStartNodes)
		updatetrucks_bns(timedelta)
		driverArcSets_online_bns()
		updateorders_bns(timedelta, currentdatetime, orderassignments)
	elseif (solutionmethod == "mvg") || (solutionmethod == "otdf") || (solutionmethod == "rrf") 
		updatedriverinfo_fragment(z_ip, currtime, driverarcstaken)
		updatedrivers_fragment(timedelta, currentdatetime, weekstart, T_off_Monday8am, currtime)
		updatetrucks(timedelta, currentdatetime, weekstart, T_off_Monday8am, x_ip, y_ip, z_ip, A_minus_irest)
		driverArcSets_online_fragment(numlocs, numarcs, numnodes, prearcs, drivers, tstep, horizon, nodes, arcs, assignedDrivers, A_minus, A_plus, T_off, drivershift, driverHomeLocs, T_off_0, shiftlength)
		yarcreduction_online(numarcs, availableDrivers, A_space, numnodes, A_plus, A_minus)
		updateorders(x_ip, y_ip, z_ip, timedelta, currentdatetime, A_minus_irest, A_plus_irest)
	elseif solutionmethod == "ptp"
		updatedriverinfo_fragment(z_ip, currtime, driverarcstaken)
		updatedrivers_fragment(timedelta, currentdatetime, weekstart, T_off_Monday8am, currtime)
		updatetrucks(timedelta, currentdatetime, weekstart, T_off_Monday8am, x_ip, y_ip, z_ip, A_minus_i_full)
		driverArcSets_online_ptp(A_minus, A_plus)
		yarcreduction_online(numarcs, availableDrivers, A_space, numnodes, A_plus, A_minus)
		updateorders_ptp(x_ip, y_ip, z_ip, timedelta, currentdatetime, orderArcSet_full, A_minus_i_full, A_plus_i_full)
	else
		updatedrivers(timedelta, currentdatetime, weekstart, T_off_Monday8am, z_ip, currtime)
		updatetrucks(timedelta, currentdatetime, weekstart, T_off_Monday8am, x_ip, y_ip, z_ip, A_minus_irest)
		driverArcSets_online(numlocs, numarcs, numnodes, prearcs, drivers, tstep, horizon, nodes, arcs, assignedDrivers, A_minus, A_plus, T_off, drivershift, driverHomeLocs, T_off_0, shiftlength)
		yarcreduction_online(numarcs, availableDrivers, A_space, numnodes, A_plus, A_minus)
		updateorders(x_ip, y_ip, z_ip, timedelta, currentdatetime, A_minus_irest, A_plus_irest)
	end

	if solutionmethod == "mvg"
		updateorderarcsets_mvg(timedelta)
		#updateorderarcsets_mvg_full(timedelta)
		updateorderarcsets_full_recalc(prearcs, shortesttriptimes)
	elseif (solutionmethod == "otd") || (solutionmethod == "otdf") || (solutionmethod == "rrf") 
		#orderarcreduction_online(prearcs, shortesttriptimes)
		updateorderarcsets_full_recalc(prearcs, shortesttriptimes)
	elseif (solutionmethod == "ptp") 
		#orderarcreduction_online(prearcs, shortesttriptimes)
		updateorderarcsets_full_recalc_ptp()
	end

	#If it's not the last iteration, get next orders
	#If it is the last iteration, print final reports
	if currtime != timedelta*(numiterations_online-1)

		getnextorders(timedelta, currentdatetime, orders, "data/lh_data_isb_connect_clean.csv", "data/vnt_data_isb_connect_clean.csv")
		
		#Order arc reduction for rivigo routes method must come after new orders are added
		if (solutionmethod == "rr") 
			orderarcreduction_online(prearcs, shortesttriptimes)
		end
			
	elseif currtime == timedelta*(numiterations_online-1)

		newordersbyiter[currtime] = 0

		#Assess delay and miles penalties for incomplete orders
		assessendofhorizonpenalties(currtime)

		#Write final output files
		if (writeresultsfile == 1) & (solutionmethod == "abcg")
			writeresults_final_online(resultsfilename, currtime)
		elseif (writeresultsfile == 1) & (solutionmethod == "otd")
			writeresults_final_online(resultsfilename, currtime)
		elseif (writeresultsfile == 1) & (solutionmethod == "rr")
			writeresults_final_online(resultsfilename, currtime)
		elseif (writeresultsfile == 1) & (solutionmethod == "bns")
			writeresults_bns(resultsfilename, currtime, "final")	
		end

		#Corridor file = number of trips between every pair of locations
		if writecorridorfile == 1
			writecorridorcounterfile(solutionmethod)
		end

		#Order outcome file = delivery time and distance for each order
		writeorderoutcomesfile(orderoutcomesfilename, solutionmethod, currtime)

		if makespatialvizfiles == 1
			spatialviz_shortestpath(string(vizfoldername, "/SpatialNetworkMaps/", vizfilename, currtime, "_shortestpaths.png"), darkestcolor, lightestcolor)
		end

		#Write file for figures = heuristic comparison
		writefiguresfile(heuristicfigurefilename)

		#COmputational summary file for abcg only
		if solutionmethod == "abcg"
			writecomputationaltimesfile(computationaltimefigurefilename, currtime)
		end

		#Save full solution
		if solutionmethod == "bns"
			savefullsolution_bns(fullxsolutionfilename, fullysolutionfilename, fullzsolutionfilename, currtime)
		end

	end
    =#

end


#---------------------------------------------------------------------------#

println("Done!")

