
using JuMP, Gurobi, Random, CSV, DataFrames, Statistics, Dates, SparseArrays 

#-----------------------------------LOAD OTHER FILES------------------------------------#

include("scripts/instancegeneration/readrivigodata.jl")
include("scripts/instancegeneration/shortestpath.jl")
include("scripts/instancegeneration/constraintmatrix.jl")
include("scripts/metrics/static_writeresults.jl")
include("scripts/arcbasedmodel/arclpandip.jl")
include("scripts/journeybasedmodel/fragmentlpandip.jl")
include("scripts/journeybasedmodel/solvedriverextensionmodel.jl")
include("scripts/journeybasedmodel/lpipbasis.jl")
include("scripts/journeybasedmodel/initializejourneymodel.jl")

include("scripts/metrics/writeresultsforearlytests.jl")

#-------------------------------------FOUR INSTANCES------------------------------------#  

inittime = time()

loclist = [20, 40, 60, 66, 66]
driverlist = [70, 300, 1600, 3316, 1600]
trucklist = [60, 230, 1200, 2495, 1200]
seedlist = [202481, 155702, 731761, 963189, 731762]
weekstartlist = [DateTime(2019, 7, 21, 8), DateTime(2019, 7, 14, 8), DateTime(2019, 7, 7, 8), DateTime(2019, 6, 30, 8), DateTime(2019, 7, 7, 8)]

ordercapslist = [3, 6, 18, 10000, 18]  #Cap on the number of orders introduced in each 6 hour increment in the online problem (6 hrs doesn't depend on tstep or timedelta to ensure consistent instance when these parameters are adjusted)
								       #Order caps is technically an online parameter, but we use it in the static version as well to ensure the orders included in the instance are the same in the static and online version

hubdistancesfilename = "data/hubdistances.csv"
traveltimesfilename = "data/traveltimes_outliers.csv"
hubdataisbfilename = "data/hub_data_isb_connect.csv"
vntdataisbfilename = "data/vnt_data_isb_connect_clean.csv"
lhdataisbfilename = "data/lh_data_isb_connect_clean.csv"

#----------------------------------INSTANCE PARAMETERS----------------------------------#  	

#Read experiment parameters 
experiment_id = ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
paramsfilename = "data/newmodel.csv"
expparms = CSV.read(paramsfilename, DataFrame)
formulation = expparms[experiment_id, 15]  # Drivers = homogeneous, heterogeneous
ex = expparms[experiment_id, 2]		
solutionmethod = expparms[experiment_id, 3]		
weekstart = expparms[experiment_id, 4]
horizon = expparms[experiment_id, 5] * 24
tstep = expparms[experiment_id, 6]
driverfactor = expparms[experiment_id, 7]
k = expparms[experiment_id, 8]   #This is the rho value
opt_gap = expparms[experiment_id, 9]
lambda = expparms[experiment_id, 10]
println("Lambda = ", lambda)

#New parameters
maxweeklydriverhours = expparms[experiment_id, 11]
lambda2 = expparms[experiment_id, 12]
variablefixingthreshold = expparms[experiment_id, 13]
strengthenedreducedcost_flag = expparms[experiment_id, 14]

#Transform date
#weekstart = DateTime(weekstart, "yyyy-mm-dd HH:MM-00")
weekstart = DateTime(weekstart) + Dates.Hour(8)

#Definition of the instance 
iterationordercap = ordercapslist[ex]					 
maxlocs = loclist[ex]
maxdrivers = round(driverlist[ex] / driverfactor, digits = 0)							 
numtrucks = round(trucklist[ex] / driverfactor, digits = 0)
randomseedval = seedlist[ex]
Random.seed!(randomseedval)

#User controlled instance parameters
#tstep = 6											# Time discretization in hours
#horizon = 24*3										# Length of planning horion in hours
shiftlength = 12									# Length of each driver shift in hours
#lambda = 500										# Objective function = lambda * delay + miles
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

#Create id for this run
runid = string("ex", ex, "_", solutionmethod, "_warmstart_exp", experiment_id, "_rundate", today())

#Travel time calculation parameters
excludeoutliers_flag = 1							# 0 = include outliers in travel time calculation, 1 = exclude outliers
googlemapstraveltimes_flag = 1						# 1 = travel time between two locations is max(Avg from Rivigo data, Google Maps travel time) (<5% of arcs use google maps time), 0 = travel time is avg of Rivigo data
includesymmetricarcs_flag = 1						# 1 = if arc A-->B present in Rivigo data but not B-->A, create synthetic arc B-->A; 0 = do not include synthetic arcs (may cause feasibility issues)
traveltimefordelay_flag = 2 						# 0 = use rounded travel times for shortest path used in delay objective, 1 = use raw travel times (best for comparing across multiple tsteps), 2 = use rounded travel times, except on the final leg of the journey where raw is used 

#MVG algorithm control parameters
newreducedcost_flag = 0								# 1 = use z-variable reduced costs to bound arc reduced costs for the subproblem, 0 = do not
fewarcsatatime_flag = 0
arcsperiteration = 5
minarcspernode = 5
globalrcthreshold = 0
consolidatepitstopsindp_flag = 0
dp_pitstopclusterradius = 70
solvedpheuristically_flag = 0
pathsperiter = 1
adaptivepathsperiter_flag = 0
pathsperiter1, pathsitercap1 = 10, 20
pathsperiter2, pathsitercap2 = 5, 20
pathsperiter3 = 1

#Uniform k and ABCG + k control parameters
#Get k directly from the solutionmethod specified above 
k = 0
ktype_flag = "pct"									# "hrs" = # of hours acceptable delay, "pct" = acceptable delay as percent of shortest path time, "min24" = max(24 hrs, percent of shortest path)

#Full IP time control parameters (also used for uniform k) 
iptimelimit = 60*95									# Max run time of IP in minutes
pbcgtimelimit = 60*95

#Branch-and-price parameters
drivervalidinequalities_flag = 0
maxbbnodes = 100000
#variableselectionstrategy_output = expparms[experiment_id, 10]
#variableselectionstrategy = split(variableselectionstrategy_output, "-") #split(expparms[experiment_id, 6],"-") 
#nodeselectionstrategy = expparms[experiment_id, 11]
#boundinitialization = expparms[experiment_id, 12]
#bap_optimality_tolerance = expparms[experiment_id, 13]

#Visualization/output/print statements control parameters
printstatements = 1									# Turn on/off printing 			
writeresultsfile = 1								# 1 = create the main results file, 0 = do not
writeorderoutcomesfiles = 1							# 1 = create the order outcomes file to see details of specific orders at specific iterations, 0 = do not
writecorridorfile = 1								# 1 = create the corridor file to see metrics on which arcs were used most often, 0 = do not
maketimespacevizfiles = 0							# Create one time-space network visualization per online iteration
makespatialvizfiles = 0							 	# Create one spatial network visualization per online iteration
makeadvancedvizfiles = 0							# Create order- and driver-specific files (LOTS of files)
darkestcolor, lightestcolor = 5, 90					# The darkest and lightest colors used on the spatial maps (0-100)
showdrivers_flag = 1								# 1 = show driver assignments on time-space network visualization, 0 = do not
saveconvergencedata_flag = 1
vizflag = 0

#File names					
vizfoldername = string("visualizations/static/run ", runid)
csvfoldername = string("outputs/bigtable_orig/")
vizfilename = string(solutionmethod)			#Folder names + file extensions added later for viz files
#resultsfilename = string(csvfoldername, runid, "/", solutionmethod, "_output.csv")
resultsfilename = string(csvfoldername, runid, "_output.csv")
orderoutcomesfilename = string(csvfoldername, runid, "/orderoutcomes.csv")
orderoutcomesbyiterfilename = string(csvfoldername, runid, "/orderoutcomes_iterationlevel.csv")
corridorfilename = string(csvfoldername, runid, "/corridorcounter.csv")
computationaltimefigurefilename = string(csvfoldername, "figures/comptimes_", runid, ".csv")
fullxsolutionfilename = string(csvfoldername, runid, "/x_soln.csv")
fullysolutionfilename = string(csvfoldername, runid, "/y_soln.csv")
fullzsolutionfilename = string(csvfoldername, runid, "/z_soln.csv")
convergencedatafilename = string(csvfoldername, "convergence_exp", runid, ".csv")
#convergencefilename = string(csvfoldername, "convergence_exp", experiment_id,".csv")


#---------------------------------IMPORT FINAL SCRIPTS----------------------------------# 

if formulation == "heterogeneous"
	include("scripts/journeybasedmodel/solvejourneymodel.jl")
	include("scripts/journeybasedmodel/solvejourneymodel_paths.jl")
	include("scripts/multiarcgeneration/multiarcgeneration.jl")
	include("scripts/columngeneration/columngeneration.jl")
elseif formulation == "homogeneous"
	include("scripts/journeybasedmodel/solvejourneymodel_homogeneous.jl")
	include("scripts/multiarcgeneration/multiarcgeneration_homogeneous.jl")
	include("scripts/columngeneration/columngeneration_homogeneous.jl")
	include("scripts/journeybasedmodel/solvejourneymodel_paths_homogeneous.jl")
end

if maketimespacevizfiles + makespatialvizfiles + makeadvancedvizfiles + vizflag >= 1
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
m_0, m_end = truckdistribution(numtrucks, numlocs, N_0, N_end)
prearcs, arcLength, arcLength_raw = readarcs(traveltimesfilename, hubdistancesfilename, tstep, numlocs, hubsTravelTimeIndex, roundup_flag, excludeoutliers_flag, hubsReverseLookup, googlemapstraveltimes_flag)
arcs, arcLookup, A_plus, A_minus, A_space, A_plus_time, A_minus_time, A_minus_space, A_plus_space, numarcs, truetraveltime = arccreation(prearcs, horizon, tstep, numnodes, nodes, numlocs)
c, u = calcobjectivecosts(hubdistancesfilename)

#Initialize orders
orderwindowstart, orderwindowend = weekstart, weekstart + Dates.Hour(horizon) - Dates.Second(1)
includeorderidlist = generateorderlist(lhdataisbfilename, vntdataisbfilename, iterationordercap, numlocs)
numorders, originloc, destloc, available, duedate, usedorderidlist, psseq, orderOriginalStartLoc, ordersinprogress = pullorders_initrivigoroutes(lhdataisbfilename, vntdataisbfilename, 10000, orderwindowstart, orderwindowend, tstep, horizon, prearcs, numlocs, tstepforordercreation, includeorderidlist)
orders = [i for i in 1:numorders]
highestorderindex = numorders
Origin, Destination = formatorders(numorders, originloc, destloc, available, duedate, tstep)
N_flow_i = flowbalancenodesets_i(orders, numnodes, Origin, Destination)

#Derive additional sets needed for various algorithms
include("scripts/instancegeneration/completeinstance.jl")
orderOriginalStartTime, orderintransit_flag = findorderstartsandtransits()
loctruckcounter, trucksintransit = findtrucksintransit()
m_0 = adjust_m_0(m_0, loctruckcounter)
driversintransit, drivers, driverStartNodes, driverEndNodes, driverHomeLocs, assignedDrivers, N_flow_t, N_flow_d, alltimeswithinview, T_off_Monday8am, T_off, drivershift, T_off_0, T_off_constr, numshifts, T_on_0 = getdriverandshiftinfo()
distbetweenlocs, shortesttriptimes, shortestpatharclists, traveltimebetweenlocs_rdd, traveltimebetweenlocs_raw, traveltimebetweenlocs_llr = findtraveltimesanddistances()
nodesLookup, arcLookup, A_minus, A_plus, c, Destination, extendednodes, extendednumnodes, extendedarcs, extendednumarcs = extendtimespacenetwork(nodesLookup, arcLookup, A_minus, A_plus, c, Destination)
arcLookup, nodesLookup, arcfinishtime, dummyarc, allarcs = calcarcfinishtimes()

#----------------------------------CREATE ARC SETS-----------------------------------# 

include("scripts/instancegeneration/initializearcsets.jl")
include("scripts/multiarcgeneration/initializeorderarcsets.jl")
primaryarcs, extendedtimearcs, orderarcs, driverarcs, hasdriverarcs = initializearcsets(A_space, A_plus, A_minus)
R_off = findreturnhomearcsets(driverarcs)
magarcs = initializeorderarcsets(k)
driversets, driverSetStartNodes, numfragments, fragmentscontaining, F_plus_ls, F_minus_ls, N_flow_ls, numeffshifts, effshift, shiftsincluded, fragdrivinghours, fragworkinghours = initializejourneymodel(maxnightsaway)

# ----------------------------------- START HERE ----------------------------------- #
# ASK YOURSELF: COULD WE DO THE SAME TRICKS WITH LP HEURISTICS THOUGH?
# ----------------------------------- START HERE ----------------------------------- #

if solutionmethod == "lp"

	lp_obj, x_lp, z_lp, lp_time, lp_bound = solvejourneymodel(1, opt_gap, orderarcs, numeffshifts)
	timeslist = (mp=lp_time, pp=0, pppar=0, ip=0)
	writeresultsforearlytests(resultsfilename, 0, "LP", lp_obj, timeslist, sum(length(orderarcs.A[i]) for i in orders))

elseif solutionmethod == "ip"

	ip_obj, x_ip, z_ip, ip_time, ip_bound = solvejourneymodel(0, opt_gap, orderarcs, numeffshifts)
	timeslist = (mp=0, pp=0, pppar=0, ip=ip_time)
	writeresultsforearlytests(resultsfilename, 0, "IP", ip_obj, timeslist, sum(length(orderarcs.A[i]) for i in orders))
	#=
	include("scripts/visualizations/timespacenetwork.jl")
	for i in orders
		arclistlist = [orderarcs.A[i], [a for a in orderarcs.A[i] if value(x_ip[i,a]) > 1e-4]]
		colorlist = [(200,200,200), (0,0,0)]
		timespacenetwork(string("outputs/viz/order", i,"_ip.png"), arclistlist, colorlist, 2000, 1200)
	end
	=#
	
	if formulation == "heterogeneous"
		println("Driver util = ", sum(sum(fragworkinghours[driverHomeLocs[d],drivershift[d],f] * value(z_ip[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers) / sum(maxweeklydriverhours for d in drivers))
	end

elseif solutionmethod == "basisip"

	lp_obj, x_lp, z_lp, lp_time, lp_bound, lpbasisarcs = solvejourneymodel(1, opt_gap, orderarcs, numeffshifts)
	bip_obj, x_bip, z_bip, bip_time, bip_bound = solvejourneymodel(0, opt_gap, lpbasisarcs, numeffshifts)

	timeslist1 = (mp=lp_time, pp=0, pppar=0, ip=0)
	writeresultsforearlytests(resultsfilename, 0, "LP", lp_obj, timeslist1, sum(length(lpbasisarcs.A[i]) for i in orders))
	timeslist2 = (mp=0, pp=0, pppar=0, ip=bip_time)
	writeresultsforearlytests(resultsfilename, 1, "IP", bip_obj, timeslist2, sum(length(lpbasisarcs.A[i]) for i in orders))
	
	#=
	include("scripts/visualizations/timespacenetwork.jl")
	for i in orders
		arclistlist = [lpbasisarcs.A[i], [a for a in lpbasisarcs.A[i] if value(x_bip[i,a]) > 1e-4]]
		colorlist = [(200,200,200), (0,0,0)]
		timespacenetwork(string("outputs/viz/order", i,"_basis.png"), arclistlist, colorlist, 2000, 1200)
	end
	=#

	if formulation == "heterogeneous"
		println("Driver util = ", sum(sum(fragworkinghours[driverHomeLocs[d],drivershift[d],f] * value(z_bip[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers) / sum(maxweeklydriverhours for d in drivers))
	end

elseif (solutionmethod == "mag") || (solutionmethod == "sag")

	magarcs = initializeorderarcsets(k)
	mag_obj, smp, x_smp, y_smp, z_smp, w_smp, magarcs, smptime, pptime, pptime_par, totalmagarcs, mag_iter = multiarcgeneration!(magarcs, variablefixingthreshold, hasdriverarcs)
	magip_obj, x_magip, z_magip, magip_time, magip_bound = solvejourneymodel(0, opt_gap, magarcs, numeffshifts)

	timeslist1 = (mp=smptime, pp=pptime, pppar=pptime_par, ip=0)
	writeresultsforearlytests(resultsfilename, 0, mag_iter, mag_obj, timeslist1, totalmagarcs)
	timeslist2 = (mp=0, pp=0, pppar=0, ip=magip_time)
	writeresultsforearlytests(resultsfilename, 1, "IP", magip_obj, timeslist2, totalmagarcs)

	#=
	include("scripts/visualizations/timespacenetwork.jl")
	for i in orders
		selectedarcs = [a for a in magarcs.A[i] if value(x_magip[i,a]) > 1e-4]
		missingarcs = [a for a in orderarcs.A[i] if (value(x_ip[i,a]) > 1e-4) & !(a in magarcs.A[i])]
		iparcs = [a for a in orderarcs.A[i] if value(x_ip[i,a]) > 1e-4]

		arclistlist = [magarcs.A[i], iparcs, missingarcs, selectedarcs]
		colorlist = [(160,160,160), (0,0,0), (0,0,240), (244,143,20)] 
		thicknesslist = [4,8,8,10]
		timespacenetwork(string("outputs/viz/order", i,"_mag.png"), arclistlist, colorlist, thicknesslist, 2400, 1800)
	end
	=#

	if formulation == "heterogeneous"
		println("Driver util = ", sum(sum(fragworkinghours[driverHomeLocs[d],drivershift[d],f] * value(z_magip[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers) / sum(maxweeklydriverhours for d in drivers))
	end

elseif solutionmethod == "cg"

	dummypath = 1
	cg_obj, rmp, x_rmp, y_rmp, z_rmp, w_rmp, cgpaths, delta, rmptime, pptime, pptime_par, totalcgpaths, cg_iter = columngeneration!(orderarcs, hasdriverarcs)
	cgip_obj, x_cgip, z_cgip, cgip_time, cgip_bound = solvejourneymodel_paths(0, opt_gap, cgpaths, delta, numeffshifts)

	timeslist1 = (mp=rmptime, pp=pptime, pppar=pptime_par, ip=0)
	writeresultsforearlytests(resultsfilename, 0, cg_iter, cg_obj, timeslist1, totalcgpaths)
	timeslist2 = (mp=0, pp=0, pppar=0, ip=cgip_time)
	writeresultsforearlytests(resultsfilename, 1, "IP", cgip_obj, timeslist2, totalcgpaths)

	if formulation == "heterogeneous"
		println("Driver util = ", sum(sum(fragworkinghours[driverHomeLocs[d],drivershift[d],f] * value(z_cgip[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers) / sum(maxweeklydriverhours for d in drivers))
	end

end

#-----------------------------------------------------------------------------#

println("Done!")
