
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
include("scripts/onlineimplementation/assessendofhorizonpenalties.jl")
include("scripts/onlineimplementation/onlinereporting.jl")
include("scripts/instancegeneration/initializearcsets.jl")
include("scripts/multiarcgeneration/initializeorderarcsets.jl")
include("scripts/visualizations/timespacenetwork.jl")
include("scripts/journeybasedmodel/solvejourneymodel_online.jl")
include("scripts/onlineimplementation/updatestate/updatepastsegments.jl")
include("scripts/onlineimplementation/updatestate/updatedrivers.jl")
include("scripts/onlineimplementation/updatestate/updatetrucks.jl")
include("scripts/onlineimplementation/updatestate/getnextorders.jl")
include("scripts/onlineimplementation/updatestate/updatearcsets.jl")
include("scripts/onlineimplementation/updatestate/updateorders.jl")
include("scripts/onlineimplementation/journeygeneration.jl")
include("scripts/visualizations/timespacenetwork.jl")

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

#----------------------------------INSTANCE PARAMETERS----------------------------------#  	

#Read experiment parameters from file
experiment_id = 1 #ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
paramsfilename = "data/relayvsptp.csv"
expparms = CSV.read(paramsfilename, DataFrame)
formulation = expparms[experiment_id, 11]  # Drivers = homogeneous, heterogeneous
ex = expparms[experiment_id, 2]		
weekstart = expparms[experiment_id, 4]
horizon = expparms[experiment_id, 5] * 24
tstep = expparms[experiment_id, 6]
driverfactor = expparms[experiment_id, 7]
k = expparms[experiment_id, 8]   #This is the rho value
opt_gap = expparms[experiment_id, 9]
lambda = expparms[experiment_id, 10]
runtype = expparms[experiment_id, 12]
operations = expparms[experiment_id, 13]
ptpvsrelay = expparms[experiment_id, 14]
println("Experiment = ", experiment_id)

#Online parameters
becomesavailablehours = expparms[experiment_id, 15]	
onlineweeks = expparms[experiment_id, 16]
timedelta = expparms[experiment_id, 17]
onlinetimehorizon = 24*7*onlineweeks		  	    			# Length of the online time horizon in hours (ex. run 1 week of online iterations), should be multiple of timedelta
numiterations_online = convert(Int64, onlinetimehorizon/timedelta)		
pruneorderarcs_flag = 1								# 1 = prune order arcs that are no longer usable after each online iteration, 0 = do not
maxorderdelayhours = onlinetimehorizon * 2

#Driver movement parameters
maxnightsaway = expparms[experiment_id, 18]
maxrepositioningdistance = expparms[experiment_id, 19]

#Read algorithm control parameters from file
solutionmethod = expparms[experiment_id, 3]		

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
finallegdistancepenalty = 2.1 #0.40						# Distance penalty assessed for orders that finish beyond the planning horizon
finallegtimepenalty = 2 #0.30							# Time/delay penalty assessed for orders that finish beyond the planning horizon
dummyendtime = 1000									# Dummy time assigned to the "beyond the horizon" nodes
driveroffdays_flag = 0
vizflag = 0
saveconvergencedata_flag = 1

#Travel time calculation parameters
excludeoutliers_flag = 1							# 0 = include outliers in travel time calculation, 1 = exclude outliers
googlemapstraveltimes_flag = ptpvsrelay==1 ? 2 : 1
includesymmetricarcs_flag = 1						# 1 = if arc A-->B present in Rivigo data but not B-->A, create synthetic arc B-->A; 0 = do not include synthetic arcs (may cause feasibility issues)
traveltimefordelay_flag = 2 						# 0 = use rounded travel times for shortest path used in delay objective, 1 = use raw travel times (best for comparing across multiple tsteps), 2 = use rounded travel times, except on the final leg of the journey where raw is used 
ensureconnectivity_flag = 1

#Uniform k and ABCG + k control parameters
ktype_flag = "pct"									# "hrs" = # of hours acceptable delay, "pct" = acceptable delay as percent of shortest path time, "min24" = max(24 hrs, percent of shortest path)

#Create id for this run
runid = string("ex", ex, "_exp", experiment_id, "_", operations, "_rundate", today())

#File names					
vizfoldername = string("visualizations/static/run ", runid)
csvfoldername = string("outputs/online/")
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
#currentdatetime = weekstart

#Create node and arc networks
timeperiods = horizon / tstep + 1
maxviztimeperiods = timeperiods 
hubCoords, hubsLookup, hubsReverseLookup, hubsTravelTimeIndex, numlocs = readlocations(hubdataisbfilename, maxlocs)
nodes, nodesLookup, N_0, N_end, numnodes = timespacentwk(numlocs, tstep, horizon)
prearcs, arcLength, arcLength_raw = readandprocessarcs(operations, traveltimesfilename, hubdistancesfilename, tstep, numlocs, hubsTravelTimeIndex, roundup_flag, excludeoutliers_flag, hubsReverseLookup, googlemapstraveltimes_flag)
arcs, arcLookup, A_plus, A_minus, A_space, A_plus_time, A_minus_time, A_minus_space, A_plus_space, numarcs, truetraveltime, arcduration, arcsbetween, arcsbetween_back = arccreation(prearcs, horizon, tstep, numnodes, nodes, numlocs)
c, u = calcobjectivecosts(hubdistancesfilename)

#Get initial state of system
currstate, includeorderidlist, drivers, driverHomeLocs, drivershift, N_flow_t, T_off_Monday8am, numshifts, originloc, destloc, orderOriginalStartLoc, orderOriginalStartTime, highestorderindex, distbetweenlocs, shortestpatharclists, traveltimebetweenlocs_rdd, traveltimebetweenlocs_raw, traveltimebetweenlocs_llr, nodesLookup, arcLookup, A_minus, A_plus, c, extendednodes, extendednumnodes, extendedarcs, extendednumarcs = getinitialstate(nodesLookup, arcLookup, A_minus, A_plus, c)

#Finish extending the time-space network
arcLookup, nodesLookup, arcfinishtime, dummyarc, allarcs = calcarcfinishtimes()
basetsn = (arcsbetween=arcsbetween, arcsbetween_back=arcsbetween_back, numlocs=numlocs, arcLookup=arcLookup, nodesLookup=nodesLookup, nodes=extendednodes, arcs=extendedarcs, numarcs=numarcs, numnodes=numnodes, horizon=horizon, tstep=tstep, extendednumarcs=extendednumarcs, extendednumnodes=extendednumnodes, A_minus=A_minus, A_plus=A_plus)
ghosttsn = createghostTSN(4)

#Create initial arc sets
include("scripts/onlineimplementation/initializecurrentstatearcs.jl")
if (operations == "ptp") & (solutionmethod == "jg") 
	currarcs, currfragments, primaryarcs, extendedtimearcs, numeffshifts, journeytime = initializecurrentstatearcs(currstate, 0);
else
	currarcs, currfragments, primaryarcs, extendedtimearcs, numeffshifts, journeytime = initializecurrentstatearcs(currstate, 1);
end
totalnumjourneys = sum(currfragments.numfragments[ds] for ds in currfragments.driversets)
println("Total journeys = ", totalnumjourneys)

total_delivtime, max_delivtime, shortestpossible_delivtime, shortestpossible_ordermiles, totalemptytrips, totalemptymiles, totaltaxitrips, totaltaximiles, ordermilespenalty, orderdelaypenalty, totalpastcost, totalordertrips, totalordermiles, ordermilesoutcomes, orderdelayoutcomes, totaldriverhours, pastordersegments, pastdriversegments_space, pastdriversegments, pastemptysegments, pasttaxisegments = initializeonlinereporting()

#---------------------------------------SOLVE----------------------------------------# 

for currtime in 0:timedelta:timedelta*(numiterations_online-1)

	println("------------------------------------BEGIN ITERATION CURRTIME = $currtime------------------------------------")

	currentdatetime = weekstart + Dates.Hour(currtime)

    #Solve current instance
    if (operations == "relay") & ((solutionmethod == "mag") || (solutionmethod == "sag"))
        #mag_obj, smp, x_smp, y_smp, z_smp, w_smp, magarcs, smptime, pptime, pptime_par, totalmagarcs = multiarcgeneration!(currstate, currfragments, currarcs)    
        ip_obj, x_ip, z_ip, w_ip, y_ip, solvetime_ip, bound_ip = solvejourneymodel(0, opt_gap, -1, currentdatetime);
		candidatejourneys, basisarcs = -1, []
	elseif (operations == "relay") & (solutionmethod == "ip") 
		ip_obj, x_ip, z_ip, w_ip, y_ip, solvetime_ip, bound_ip = solvejourneymodel(0, opt_gap, -1, currentdatetime);
		candidatejourneys, basisarcs = -1, []
	elseif (operations == "relay") & (solutionmethod == "basisip") 
		lp_obj, x_lp, z_lp, w_lp, y_lp, solvetime_lp, bound_lp, basisarcs = solvejourneymodel_relayred(1, opt_gap, -1, currentdatetime, currarcs.orderarcs);
		ip_obj, x_ip, z_ip, w_ip, y_ip, solvetime_ip, bound_ip = solvejourneymodel_relayred(0, opt_gap, -1, currentdatetime, basisarcs);
		candidatejourneys = -1
	elseif (operations == "ptp") & (solutionmethod == "basisip") 
		lp_obj, x_lp, z_lp, w_lp, y_lp, solvetime_lp, bound_lp, candidatejourneys = solvejourneymodel(1, opt_gap, -1, currentdatetime);
		ip_obj, x_ip, z_ip, w_ip, y_ip, solvetime_ip, bound_ip = solvejourneymodel(0, opt_gap, candidatejourneys, currentdatetime);
		basisarcs = []
	elseif (operations == "ptp") & (solutionmethod == "ip") 
		ip_obj, x_ip, z_ip, w_ip, y_ip, solvetime_ip, bound_ip = solvejourneymodel(0, opt_gap, -1, currentdatetime);
		candidatejourneys, basisarcs = -1, []
	elseif (operations == "ptp") & (solutionmethod == "jg") 
		include("scripts/onlineimplementation/journeygeneration.jl")
		lp_obj, x_lp, z_lp, w_lp, y_lp = journeygeneration(opt_gap, currarcs.orderarcs, currarcs.ghostdriverarcs)
		candidatejourneys, basisarcs = -1, []
		ip_obj, x_ip, z_ip, w_ip, y_ip, solvetime_ip, bound_ip = solvejourneymodel(0, opt_gap, candidatejourneys, currentdatetime);
	end

	#Visualize
	#=for i in currstate.orders
		availarcs, usedarcs = [a for a in currarcs.orderarcs.A[i]], [a for a in currarcs.orderarcs.A[i] if x_ip[i,a]>1e-4]
		timespacenetwork(string("outputs/viz/online/iter",currtime,"_order",i,".png"), [availarcs, usedarcs], [(150,150,150),(0,0,0)], [3,8], ["solid", "solid"], [0,0], 2400, 1800)
	end=#

	#Find arcs taken by drivers
	driverarcstaken = updatepastsegments(timedelta, x_ip, y_ip, z_ip, w_ip, candidatejourneys, currentdatetime, basisarcs)
	updatelasttimehome(driverarcstaken)
	
	#Update local datetime 
	currentdatetime = currentdatetime + Dates.Hour(timedelta)

	#Iterate forward by timedelta 
	updatedriverlocations(currentdatetime, driverarcstaken);
	updatedriversshifts(currentdatetime, weekstart, T_off_Monday8am);
	updatetrucks(timedelta, currentdatetime, weekstart, x_ip, y_ip, basisarcs);
	updatedriverarcsets();

	#=myarcs = [a for a in currarcs.driverarcs.A[driverHomeLocs[7],drivershift[7]]]
	include("scripts/visualizations/timespacenetwork.jl")
	timespacenetwork("outputs/viz/aaa_all.png", [myarcs], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)

	myarcs=[]
	for j in newfragments.fragments[17, 1, 17, 0]
		myarcs = union(myarcs,newfragments.fragmentarcs[17, 1, 17, 0, j])
	end
	timespacenetwork("outputs/viz/aaa_all.png", [myarcs], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)=#

	#Better way to do this without the global var update?
	newfragments = updatedriverjourneys(1)
	global currfragments = newfragments

	#Update orders
	updateorders(x_ip, timedelta, currentdatetime, basisarcs)
	
	#Better way to do this without the global var update?
	neworderarcs, newmagarcs = updateorderarcs()
	global currarcs = (orderarcs=neworderarcs, driverarcs=currarcs.driverarcs, hasdriverarcs=currarcs.hasdriverarcs, magarcs=newmagarcs, ghostdriverarcs=currarcs.ghostdriverarcs) 

	#If it's not the last iteration, get next orders
	if currtime != timedelta*(numiterations_online-1)

		getnextorders(timedelta, currentdatetime, lhdataisbfilename, vntdataisbfilename)
		writeresults_onlineiteration(resultsfilename, currtime)

	#If it is the last iteration, print final reports
	elseif currtime == timedelta*(numiterations_online-1)

		writeresults_onlineiteration(resultsfilename, currtime)

		#Assess delay and miles penalties for incomplete orders
		assessendofhorizonpenalties(currstate, currtime)

		writeresults_onlinefinal(resultsfilename, currtime + tstep)

	end

end

println("Total past cost for exp. $experiment_id = ", sum(values(totalpastcost)))

#---------------------------------------------------------------------------#

#= 
for i in currstate.orders
	myarcs = [a for a in currarcs.orderarcs.A[i]]
	timespacenetwork(string("outputs/viz/aaa_order",i,".png"), [myarcs], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
end
for item in currfragments.driversets
	myarcs = []
	for f in currfragments.fragments[item]
		myarcs = union(myarcs, currfragments.fragmentarcs[item[1], item[2], item[3], item[4], f])
	end
	timespacenetwork(string("outputs/viz/aaa_driver",item,".png"), [myarcs], [(150,150,150)], [3], ["solid"], [0], 2400, 1800)
end
=#

#---------------------------------------------------------------------------#

println("Done!")
