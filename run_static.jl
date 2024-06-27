
using JuMP, Gurobi, Random, CSV, DataFrames, Statistics, Dates, SparseArrays 

#-----------------------------------LOAD OTHER FILES------------------------------------#

include("scripts/instancegeneration/readrivigodata.jl")
include("scripts/instancegeneration/shortestpath.jl")
include("scripts/instancegeneration/constraintmatrix.jl")
include("scripts/instancegeneration/completeinstance.jl")
include("scripts/instancegeneration/initializearcsets.jl")
include("scripts/multiarcgeneration/initializeorderarcsets.jl")
include("scripts/journeybasedmodel/initializejourneymodel.jl")
include("scripts/metrics/writeresultsforrun.jl")
include("scripts/onlineimplementation/initializecurrentstatearcs.jl")

#-------------------------------------FOUR INSTANCES------------------------------------#  

inittime = time()
starttime = time()

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
experiment_id = ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
paramsfilename = "data/laborandshift_sensitivity.csv"
expparms = CSV.read(paramsfilename, DataFrame)
formulation = expparms[experiment_id, 15]  # Drivers = homogeneous, heterogeneous
ex = expparms[experiment_id, 2]		
weekstart = expparms[experiment_id, 4]
horizon = expparms[experiment_id, 5] * 24
tstep = expparms[experiment_id, 6]
driverfactor = expparms[experiment_id, 7]
k = expparms[experiment_id, 8]  		   # This is the rho value
opt_gap = expparms[experiment_id, 9]
lambda = expparms[experiment_id, 10]
maxweeklydriverhours = expparms[experiment_id, 11]
lambda2 = expparms[experiment_id, 12]
runtype = "static"
operations = "relay" 					   # "ptp" or "relay"
ptpvsrelay = 1
println("Experiment = ", experiment_id)

#Manual parameters for response/appendix experiments
roundeddrivinghours_flag = 0
if formulation == "heterogeneous"
	csvfoldername = string("outputs/laborandshiftsensitivity/")
	deadlines_flag = 1
	finallegdistancepenalty = 0.80				# Distance penalty assessed for orders that finish beyond the planning horizon
	finallegtimepenalty = 0.70					# Time/delay penalty assessed for orders that finish beyond the planning horizon
	deadlineasmultipleofshortestpath = 2
elseif formulation == "homogeneous"
	csvfoldername = string("outputs/driverhiring/")
	deadlines_flag = 0
	finallegdistancepenalty = 0.40 			    # Distance penalty assessed for orders that finish beyond the planning horizon
	finallegtimepenalty = 0.30					# Time/delay penalty assessed for orders that finish beyond the planning horizon
	deadlineasmultipleofshortestpath = 2
elseif formulation == "homogeneousdeadlines"
	#csvfoldername = string("outputs/ordersensitivity/")
	csvfoldername = string("outputs/laborandshiftsensitivity/")
	deadlines_flag = 1
	finallegdistancepenalty = 0.40 			    # Distance penalty assessed for orders that finish beyond the planning horizon
	finallegtimepenalty = 0.30	
	deadlineasmultipleofshortestpath = expparms[experiment_id, 24]
else
	throw(DomainError(formulation, "formulation not recognized: must be 'homogeneous' or 'heterogenous'"))
end	

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
mip_focus = expparms[experiment_id, 23]  
driverstohire = expparms[experiment_id, 24]  

#Transform date
weekstart = DateTime(weekstart) + Dates.Hour(8)

#Definition of the instance 
iterationordercap = ordercapslist[ex]					 
maxlocs = loclist[ex]
maxdrivers1 = round(driverlist[ex] / driverfactor, digits = 0)	
maxdrivers2 = round(driverlist[ex] / driverfactor, digits = 0)	
hiredriversto = "all"						 
numtrucks = round(trucklist[ex] / driverfactor, digits = 0)
randomseedval = seedlist[ex]
Random.seed!(randomseedval)

#Problem/optimization parameters (all hard coded values we decided on at the beginning)
shiftlength = 12									# Length of each driver shift in hours
taxicostpct = 2.0                          			# Cost to taxi along an arc as a percentage of cost to drive a truck along that arc (should be > 1.0)
roundup_flag = 1			 					 	# 0 = round down to find discretized travel times, 1 = round up
drivershifttstep = 12								# How many hours between start of driver shifts, (ex. drivershifttstep=12 means each driver's first shift starts at time 0 or time 12, drivershifttstep=6 means a driver's first shift could start at time 0, 6, 12, or 18)
tstepforordercreation = 12 							# Should be same as timedelta in equivalent online instance; used to round order available time stamps from Rivigo data (ex. round observed available time to previous 12 hour block)
inprogressdummyarc_flag = 0					    	# 1 = allow in progress orders to be assigned to the dummy arc, 0 = do not (Should be assigned to 0 to ensure feasibility/continuity of online iterations)
truearcfinishtime_flag = 0							# 1 = use unrounded arc travel times to assess order delivery delay (still has some bugs), 0 = use travel times rounded up to the next time step
dummyendtime = 1000									# Dummy time assigned to the "beyond the horizon" nodes
maxnightsaway = 1
driveroffdays_flag = 0
vizflag = 0
saveconvergencedata_flag = 1
timedeltaexp_flag = 0
ordergenerationtstep = 48

#Travel time calculation parameters
excludeoutliers_flag = 1							# 0 = include outliers in travel time calculation, 1 = exclude outliers
googlemapstraveltimes_flag = ptpvsrelay==1 ? 2 : 1	# 2 = travel time between two locations is Google maps, 1 = max(Avg from Rivigo data, Google Maps travel time) (<5% of arcs use google maps time), 0 = travel time is avg of Rivigo data
includesymmetricarcs_flag = 1						# 1 = if arc A-->B present in Rivigo data but not B-->A, create synthetic arc B-->A; 0 = do not include synthetic arcs (may cause feasibility issues)
traveltimefordelay_flag = 2 						# 0 = use rounded travel times for shortest path used in delay objective, 1 = use raw travel times (best for comparing across multiple tsteps), 2 = use rounded travel times, except on the final leg of the journey where raw is used 
ensureconnectivity_flag = 1
onlinetimehorizon = horizon*2

#Sensitivity analysis
percentnightshift = paramsfilename == "data/laborandshift_sensitivity.csv" ? expparms[experiment_id, 25] : 0.50
laborcost_delta = paramsfilename == "data/laborandshift_sensitivity.csv" ? expparms[experiment_id, 26] : 0
driverinventorycost_theta = paramsfilename == "data/laborandshift_sensitivity.csv" ? expparms[experiment_id, 27] : 0
basedriverfactor = driverfactor 					#Online use only 
hiredriversto = "all"								#Online use only 
timedeltaexp_flag = 0								#Online use only 

#Uniform k and ABCG + k control parameters
ktype_flag = "pct"									# "hrs" = # of hours acceptable delay, "pct" = acceptable delay as percent of shortest path time, "min24" = max(24 hrs, percent of shortest path)

#Create id for this run
runid = string("ex", ex, "_exp", experiment_id, "_", solutionmethod, "_rundate", today())

#File names					
vizfoldername = string("visualizations/static/run ", runid)
resultsfilename = string(csvfoldername, runid, "_output.csv")
convergencedatafilename = string(csvfoldername, "convergence_exp", runid, ".csv")

#---------------------------------IMPORT FINAL SCRIPTS----------------------------------# 

#Import algorithm functions
if formulation == "heterogeneous"
	include("scripts/journeybasedmodel/solvejourneymodel.jl")
	include("scripts/journeybasedmodel/solvejourneymodel_paths.jl")
	include("scripts/multiarcgeneration/multiarcgeneration_heterogeneous.jl")
	include("scripts/columngeneration/columngeneration.jl")
	include("scripts/arcbasedmodel/solvearcbasedmodel_heterogeneous.jl")
elseif (formulation == "homogeneous") || (formulation == "homogeneousdeadlines")
	include("scripts/journeybasedmodel/solvejourneymodel_homogeneous.jl")
	include("scripts/multiarcgeneration/multiarcgeneration_homogeneous.jl")
	include("scripts/columngeneration/columngeneration_homogeneous.jl")
	include("scripts/journeybasedmodel/solvejourneymodel_paths_homogeneous.jl")
	include("scripts/arcbasedmodel/solvearcbasedmodel_homogeneous.jl")
end

#Import cut functions
include("scripts/knapsackcuts/findknapsackcuts.jl")
include("scripts/knapsackcuts/solvelpwithcuts.jl")
include("scripts/knapsackcuts/solveipwithcuts.jl")

#Import extensions
include("scripts/extensions/calcorderdeadlines.jl")

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
m_0, m_end = truckdistribution(numtrucks, numlocs, N_0, N_end)
prearcs, arcLength, arcLength_raw = readandprocessarcs(operations, traveltimesfilename, hubdistancesfilename, tstep, numlocs, hubsTravelTimeIndex, roundup_flag, excludeoutliers_flag, hubsReverseLookup, googlemapstraveltimes_flag)
arcs, arcLookup, A_plus, A_minus, A_space, A_plus_time, A_minus_time, A_minus_space, A_plus_space, numarcs, truetraveltime, arcduration, arcsbetween, arcsbetween_back = arccreation(prearcs, horizon, tstep, numnodes, nodes, numlocs)
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
orderOriginalStartTime, orderintransit_flag = findorderstartsandtransits(orders, Origin)
loctruckcounter, trucksintransit = findtrucksintransit(ordersinprogress, originloc, available)
m_0 = adjust_m_0(m_0, loctruckcounter)
driversintransit, drivers, driverStartNodes, driverEndNodes, driverHomeLocs, assignedDrivers, N_flow_t, N_flow_d, alltimeswithinview, T_off_Monday8am, T_off, drivershift, T_off_0, T_off_constr, numshifts, T_on_0 = getdriverandshiftinfo()
distbetweenlocs, shortesttriptimes, shortestpatharclists, traveltimebetweenlocs_rdd, traveltimebetweenlocs_raw, traveltimebetweenlocs_llr = findtraveltimesanddistances(orders, Origin, Destination)
#println(weekstart)
#println(sum(distbetweenlocs[originloc[i], destloc[i]] for i in orders))
orderdeadline = calcorderdeadlines(shortesttriptimes)
nodesLookup, arcLookup, A_minus, A_plus, c, extendednodes, extendednumnodes, extendedarcs, extendednumarcs, u = extendtimespacenetwork(nodesLookup, arcLookup, A_minus, A_plus, c, u, distbetweenlocs)
Destination = extendDestination(orders, Destination, extendednodes)
arcLookup, nodesLookup, arcfinishtime, dummyarc, allarcs = calcarcfinishtimes()

#----------------------------------CREATE ARC SETS-----------------------------------# 

basetsn = (arcsbetween=arcsbetween, arcsbetween_back=arcsbetween_back, numlocs=numlocs, arcLookup=arcLookup, nodesLookup=nodesLookup, nodes=extendednodes, arcs=extendedarcs, numarcs=numarcs, numnodes=numnodes, horizon=horizon, tstep=tstep, extendednumarcs=extendednumarcs, extendednumnodes=extendednumnodes, A_minus=A_minus, A_plus=A_plus)
ghosttsn = createghostTSN(maxnightsaway+2)
lasttimehome = [0 for d in drivers]

currstate = (m_0=m_0, m_end=m_end, trucksintransit=trucksintransit, orders=orders, 
    Origin=Origin, Destination=Destination, driversintransit=driversintransit, N_flow_d=N_flow_d, N_flow_i=N_flow_i,
    alltimeswithinview=alltimeswithinview, T_off=T_off, T_off_0=T_off_0, T_off_constr=T_off_constr, T_on_0=T_on_0,
    available=available, duedate=duedate, usedorderidlist=usedorderidlist, psseq=psseq, ordersinprogress=ordersinprogress, 
    shortesttriptimes=shortesttriptimes, orderintransit_flag=orderintransit_flag,
    driverStartNodes=driverStartNodes, driverEndNodes=driverEndNodes, assignedDrivers=assignedDrivers,
    lasttimehome=lasttimehome)

primaryarcs, extendedtimearcs, orderarcs, driverarcs, hasdriverarcs, ghostdriverarcs = initializearcsets(A_space, A_plus, A_minus, orders, Origin, Destination, driverStartNodes, T_off)
R_off = findreturnhomearcsets(driverarcs, T_off_constr)
magarcs = initializeorderarcsets(k, orders, originloc, destloc, Origin, Destination, shortesttriptimes)
driversets, driverSetStartNodes, numfragments, fragmentscontaining, F_plus_ls, F_minus_ls, N_flow_ls, numeffshifts, effshift, shiftsincluded, fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway, fragmentarcs = initializejourneymodel(maxnightsaway, T_off, T_on_0)
nocuts=(vars=[], rhs=[], coeff=[])

#Sensitivity: labor costs and inventory costs
if paramsfilename == "data/laborandshift_sensitivity.csv"
	include("scripts/extensions/calculatelaborandinventorycosts.jl")
	journeylaborcost, journeyinventorycost = calculatelaborandinventorycosts()
end

#include("scripts/onlineimplementation/initializecurrentstatearcs.jl")
#driversets, driversingroup, numdrivergroups, drivergroupnum, drivergroupdesc, numeffshifts, effshift, shiftsincluded = finddriversets_online(T_off, driverStartNodes, lasttimehome) 
#basetsn = (arcsbetween=arcsbetween, arcsbetween_back=arcsbetween_back, numlocs=numlocs, arcLookup=arcLookup, nodesLookup=nodesLookup, nodes=extendednodes, arcs=extendedarcs, numarcs=numarcs, numnodes=numnodes, horizon=horizon, tstep=tstep, extendednumarcs=extendednumarcs, extendednumnodes=extendednumnodes, A_minus=A_minus, A_plus=A_plus)
#ghosttsn = createghostTSN(maxnightsaway)

#journeystart = time()
#numfragments, fragmentscontaining, fragmentarcs, F_plus_g, F_minus_g, N_flow_g = initializedriversetjourneys(driversets, drivergroupnum, driversingroup, drivergroupdesc, driverarcs, ghostdriverarcs, 1)
#journeytime = time() - journeystart
#fragdrivinghours, fragworkinghours, workingfragments, fragmentnightsaway = getfragmentstats(driversets, numfragments, fragmentarcs)

#currarcs = (orderarcs=orderarcs, driverarcs=driverarcs, hasdriverarcs=hasdriverarcs, magarcs=magarcs); #, R_off=R_off)
#currfragments = (driversets=driversets, driversingroup=driversingroup, numdrivergroups=numdrivergroups, drivergroupnum=drivergroupnum, drivergroupdesc=drivergroupdesc, numeffshifts=numeffshifts, effshift=effshift, shiftsincluded=shiftsincluded,numfragments=numfragments, fragmentscontaining=fragmentscontaining, fragmentarcs=fragmentarcs, F_plus_g=F_plus_g, F_minus_g=F_minus_g, N_flow_g=N_flow_g, fragdrivinghours=fragdrivinghours, fragworkinghours=fragworkinghours, workingfragments=workingfragments, fragmentnightsaway=fragmentnightsaway);

#println("Nights away = ", maxnightsaway)
#println("Journeys = ", sum(sum(numfragments[l,s] for s in 1:numeffshifts) for l in 1:numlocs))

#---------------------------------------SOLVE----------------------------------------# 

if solutionmethod == "lp"

	lp_obj, x_lp, z_lp, lp_time, lp_bound = solvejourneymodel(1, opt_gap, orderarcs, numeffshifts, nocuts)
	timeslist = (mp=lp_time, pp=0, pppar=0, ip=0, cut=0, full=0)
	writeresultsforrun(resultsfilename, 0, "LP", lp_obj, timeslist, sum(length(orderarcs.A[i]) for i in orders), x_lp, z_lp)

	#=include("scripts/visualizations/timespacenetwork.jl")
	for i in orders
		lparcs = [a for a in orderarcs.A[i] if value(x_lp[i,a]) > 1e-4]

		fractlist = [[],Dict(),[]]
		for a in lparcs
			fractlist[2][a] = string(round(value(x_lp[i,a]), digits=2))
		end

		arclistlist = [orderarcs.A[i], lparcs]
		colorlist = [(190,190,190), (0,0,0)] 
		thicknesslist = [3,11]
		timespacenetwork(string("outputs/viz/order", i,"_lp.png"), arclistlist, colorlist, thicknesslist, fractlist, 2400, 1800)
	end=#
	
elseif solutionmethod == "lpcuts"

	lpc_obj, x_lpc, z_lpc, lpc_time, lpc_bound, knapsackcuts, basisarcs = solvelpwithcuts(opt_gap, orderarcs, knapsackcuttype)
	timeslist = (mp=lpc_time, pp=0, pppar=0, ip=0, cut=0, full=0)
	writeresultsforrun(resultsfilename, 0, "LP", lpc_obj, timeslist, sum(length(orderarcs.A[i]) for i in orders), x_lpc, z_lpc)

	ipc_obj, x_ipc, z_ipc, ipc_time, ipc_bound = solvejourneymodel(0, opt_gap, orderarcs, numeffshifts, knapsackcuts)
	timeslist = (mp=0, pp=0, pppar=0, ip=ipc_time, cut=0, full=0)
	writeresultsforrun(resultsfilename, 1, "IP", ipc_obj, timeslist, sum(length(orderarcs.A[i]) for i in orders), x_ipc, z_ipc)

	#=include("scripts/visualizations/timespacenetwork.jl")
	for i in orders
		lparcs = [a for a in orderarcs.A[i] if value(x_lpc[i,a]) > 1e-4]
		iparcs = [a for a in orderarcs.A[i] if value(x_ip[i,a]) > 1e-4]
		arclistlist = [orderarcs.A[i], lparcs, iparcs]
		colorlist = [(200,200,200),(0,0,0), (255,0,0)] 
		thicknesslist = [5,11,6]
		fractlist = [[],Dict(),[]]
		for a in lparcs
			fractlist[2][a] = string(round(value(x_lpc[i,a]), digits=2))
		end
		timespacenetwork(string("outputs/viz/order", i,".png"), arclistlist, colorlist, thicknesslist, fractlist, 2000, 1200)
	end
	for d in drivers
		lparcs, iparcs = [], []
		fractlist = [[],Dict(),[]]
		for f in 1:numfragments[driverHomeLocs[d], drivershift[d]] 
			if value(z_lpc[d,f]) > 1e-4
				for a in 1:numarcs
					if f in fragmentscontaining[driverHomeLocs[d],drivershift[d],a]
						push!(lparcs, a)
						try
							fractlist[2][a] += round(value(z_lpc[d,f]), digits=2)
						catch
							fractlist[2][a] = round(value(z_lpc[d,f]), digits=2)
						end
					end
				end
			end
		end	
		for f in 1:numfragments[driverHomeLocs[d], drivershift[d]]
			if value(z_ip[d,f]) > 1e-4
				for a in 1:numarcs
					if f in fragmentscontaining[driverHomeLocs[d],drivershift[d],a]
						push!(iparcs, a)
					end
				end
			end
		end		
		for a in keys(fractlist[2])
			fractlist[2][a] = string(round(fractlist[2][a], digits=2))
		end
		arclistlist = [driverarcs.A[d], lparcs, iparcs]
		colorlist = [(200,200,200),(0,0,0), (255,0,0)] 
		thicknesslist = [5,11,6]
		timespacenetwork(string("outputs/viz/driver", d,".png"), arclistlist, colorlist, thicknesslist, fractlist, 2000, 1200)
	end=#

elseif (solutionmethod == "ip") #& (knapsackcuttype == 0)

	lp_obj, x_lp, z_lp, lp_time, lp_bound, lpbasisarcs = solvejourneymodel(1, opt_gap, orderarcs, numeffshifts, nocuts)
	timeslist = (mp=lp_time, pp=0, pppar=0, ip=0, cut=0, full=0)
	writeresultsforrun(resultsfilename, 0, "LP", lp_obj, timeslist, sum(length(orderarcs.A[i]) for i in orders), x_lp, z_lp)

	ip_obj, x_ip, z_ip, ip_time, ip_bound = solvejourneymodel(0, opt_gap, orderarcs, numeffshifts, nocuts)
	timeslist = (mp=0, pp=0, pppar=0, ip=ip_time, cut=0, full=0)
	writeresultsforrun(resultsfilename, 1, "IP", ip_obj, timeslist, sum(length(orderarcs.A[i]) for i in orders), x_ip, z_ip)
	
	#=
	include("scripts/visualizations/timespacenetwork.jl")
	for i in orders
		arclistlist = [orderarcs.A[i], [a for a in orderarcs.A[i] if value(x_ip[i,a]) > 1e-4]]
		colorlist = [(200,200,200), (0,0,0)]
		timespacenetwork(string("outputs/viz/order", i,"_ip.png"), arclistlist, colorlist, 2000, 1200)
	end
	=#

elseif (solutionmethod == "ip") & (knapsackcuttype != 0)

	#Solve IP with generated cuts
	ipc_obj, x_ipc, z_ipc, ipc_time, ipc_bound = solveipwithcuts(opt_gap, orderarcs, numeffshifts)
	timeslist = (mp=0, pp=0, pppar=0, ip=ipc_time, cut=0, full=0)
	writeresultsforrun(resultsfilename, 1, "IP", ipc_obj, timeslist, sum(length(orderarcs.A[i]) for i in orders), x_ipc, z_ipc)

	#=
	include("scripts/visualizations/timespacenetwork.jl")
	for i in orders
		arclistlist = [orderarcs.A[i], [a for a in orderarcs.A[i] if value(x_ip[i,a]) > 1e-4]]
		colorlist = [(200,200,200), (0,0,0)]
		timespacenetwork(string("outputs/viz/order", i,"_ip.png"), arclistlist, colorlist, 2000, 1200)
	end
	=#

elseif solutionmethod == "ipred"

	reducedarcs = initializeorderarcsets(k, orders, originloc, destloc, Origin, Destination, shortesttriptimes)
	ipk_obj, x_ipk, z_ipk, ipk_time, ipk_bound = solvejourneymodel(0, opt_gap, reducedarcs, numeffshifts, nocuts)
	timeslist = (mp=0, pp=0, pppar=0, ip=ipk_time, cut=0, full=0)
	writeresultsforrun(resultsfilename, 0, "IPreduced", ipk_obj, timeslist, sum(length(reducedarcs.A[i]) for i in orders), x_ipk, z_ipk)
		
elseif solutionmethod == "arcip"

	arcip_obj, x_arcip, z_arcip, arcip_time, arcip_bound = solvearcbasedmodel(orderarcs, 0)
	timeslist = (mp=0, pp=0, pppar=0, ip=arcip_time, cut=0, full=0)
	writeresultsforrun(resultsfilename, 0, "ArcIP", arcip_obj, timeslist, sum(length(orderarcs.A[i]) for i in orders), x_arcip, z_arcip)
		
elseif (solutionmethod == "basisip") & (formulation == "homogeneous")

	lp_obj, x_lp, z_lp, lp_time, lp_bound, lpbasisarcs = solvejourneymodel(1, opt_gap, orderarcs, numeffshifts, nocuts)
	bip_obj, x_bip, z_bip, bip_time, bip_bound = solvejourneymodel(0, opt_gap, lpbasisarcs, numeffshifts, nocuts)

	timeslist1 = (mp=lp_time, pp=0, pppar=0, ip=0, cut=0, full=0)
	writeresultsforrun(resultsfilename, 0, "LP", lp_obj, timeslist1, sum(length(lpbasisarcs.A[i]) for i in orders), x_lp, z_lp)
	timeslist2 = (mp=0, pp=0, pppar=0, ip=bip_time, cut=0, full=0)
	writeresultsforrun(resultsfilename, 1, "IP", bip_obj, timeslist2, sum(length(lpbasisarcs.A[i]) for i in orders), x_bip, z_bip)
	
	println("IP-LP gap = ", round(100*(bip_obj-lp_obj)/lp_obj, digits=2), "%")

elseif solutionmethod == "basisip"

	lp_obj, x_lp, z_lp, lp_time, lp_bound, knapsackcuts, lpbasisarcs = solvelpwithcuts(opt_gap, orderarcs, knapsackcuttype)
	bip_obj, x_bip, z_bip, bip_time, bip_bound = solvejourneymodel(0, opt_gap, lpbasisarcs, numeffshifts, knapsackcuts)

	timeslist1 = (mp=lp_time, pp=0, pppar=0, ip=0, cut=0, full=0)
	writeresultsforrun(resultsfilename, 0, "LP", lp_obj, timeslist1, sum(length(lpbasisarcs.A[i]) for i in orders), x_lp, z_lp)
	timeslist2 = (mp=0, pp=0, pppar=0, ip=bip_time, cut=0, full=0)
	writeresultsforrun(resultsfilename, 1, "IP", bip_obj, timeslist2, sum(length(lpbasisarcs.A[i]) for i in orders), x_bip, z_bip)
	
	println("IP-LP gap = ", round(100*(bip_obj-lp_obj)/lp_obj, digits=2), "%")

	#=
	include("scripts/visualizations/timespacenetwork.jl")
	for i in orders
		arclistlist = [lpbasisarcs.A[i], [a for a in lpbasisarcs.A[i] if value(x_bip[i,a]) > 1e-4]]
		colorlist = [(200,200,200), (0,0,0)]
		timespacenetwork(string("outputs/viz/order", i,"_basis.png"), arclistlist, colorlist, 2000, 1200)
	end
	=#

elseif ((solutionmethod == "mag") || (solutionmethod == "sag")) & ((formulation == "homogeneous") || (formulation == "homogeneousdeadlines"))

	mag_obj, smp, x_smp, y_smp, z_smp, w_smp, magarcs, smptime, pptime, pptime_par, totalmagarcs, mag_iter = multiarcgeneration!(magarcs, hasdriverarcs) 
	magip_obj, x_magip, z_magip, magip_time, magip_bound, ipbasisarcs, totalmiles, totaldelay, totalordertime, totalemptymiles, totalrepomiles, totalshortestpathmiles, totalordermiles, totalpenaltymiles, timelaboring, timeasinventory = solvejourneymodel(0, opt_gap, magarcs, numeffshifts, nocuts)

	timeslist1 = (mp=smptime, pp=pptime, pppar=pptime_par, ip=0, cut=0, full=0)
	writeresultsforrun_deadlines(resultsfilename, 0, mag_iter, mag_obj, timeslist1, totalmagarcs, x_smp, z_smp, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
	timeslist2 = (mp=0, pp=0, pppar=0, ip=magip_time, cut=0, full=0)
	writeresultsforrun_deadlines(resultsfilename, 1, "IP", magip_obj, timeslist2, totalmagarcs, x_magip, z_magip, totalmiles, totaldelay, totalordertime, totalemptymiles, totalrepomiles, totalshortestpathmiles, totalordermiles, totalpenaltymiles, timelaboring, timeasinventory)
	
elseif (solutionmethod == "mag") || (solutionmethod == "sag")

	variableusecount, startercuts, starterfixedvars = initmagsets(magarcs)	
	mag_obj, smp, x_smp, y_smp, z_smp, w_smp, magarcs, smptime, pptime, pptime_par, totalmagarcs, mag_iter, knapsackcuts, cuttime, arcsbeforecolmgmt, arcsbeforevarsetting = multiarcgeneration_heterogeneous!(magarcs, hasdriverarcs, startercuts, starterfixedvars, variableusecount, 0, 1, 100000)
	
	#mag_obj, smp, x_smp, y_smp, z_smp, w_smp, magarcs, smptime, pptime, pptime_par, totalmagarcs, mag_iter, knapsackcuts, cuttime = multiarcgeneration!(magarcs, variablefixingthreshold, hasdriverarcs)

	magip_obj, x_magip, z_magip, magip_time, magip_bound = solvejourneymodel(0, opt_gap, magarcs, numeffshifts, nocuts)

	timeslist1 = (mp=smptime, pp=pptime, pppar=pptime_par, ip=0, cut=cuttime, full=0)
	writeresultsforrun(resultsfilename, 0, mag_iter, mag_obj, timeslist1, totalmagarcs, x_smp, z_smp)
	timeslist2 = (mp=0, pp=0, pppar=0, ip=magip_time, cut=0, full=0)
	writeresultsforrun(resultsfilename, 1, "IP", magip_obj, timeslist2, totalmagarcs, x_magip, z_magip)
	
	writedriverstats(string("outputs/bigtable_new/driverstats_exp", experiment_id,".csv"), z_magip)
	
	#=
	include("scripts/visualizations/timespacenetwork.jl")
	for i in orders
		magiparcs = [a for a in magarcs.A[i] if value(x_magip[i,a]) > 1e-4]
		#iparcs = [a for a in orderarcs.A[i] if value(x_ip[i,a]) > 1e-4]
		#missingarcs = [a for a in orderarcs.A[i] if (value(x_ip[i,a]) > 1e-4) & !(a in magarcs.A[i])]
		lparcs = [a for a in orderarcs.A[i] if value(x_lp[i,a]) > 1e-4]

		fractlist = [[],Dict(),[]]
		for a in lparcs
			fractlist[2][a] = string(round(value(x_smp[i,a]), digits=2))
		end

		arclistlist = [magarcs.A[i], lparcs, magiparcs]
		colorlist = [(190,190,190), (0,0,0), (254,97,0)] 
		thicknesslist = [3,11,7]
		timespacenetwork(string("outputs/viz/order", i,"_mag.png"), arclistlist, colorlist, thicknesslist, fractlist, 2400, 1800)
	end
	for d in drivers
		lparcs, iparcs = [], []
		for f in 1:numfragments[driverHomeLocs[d], drivershift[d]]
			if value(z_magip[d,f]) > 1e-4
				for a in 1:numarcs
					if f in fragmentscontaining[driverHomeLocs[d],drivershift[d],a]
						push!(iparcs, a)
					end
				end
			end
		end
		for f in 1:numfragments[driverHomeLocs[d], drivershift[d]]
			if value(z_smp[d,f]) > 1e-4
				for a in 1:numarcs
					if f in fragmentscontaining[driverHomeLocs[d],drivershift[d],a]
						push!(lparcs, a)
					end
				end
			end
		end
		
		arclistlist = [driverarcs.A[d], lparcs, iparcs]
		colorlist = [(200,200,200),(0,0,0), (254,97,0)] 
		thicknesslist = [3,11,6]
		timespacenetwork(string("outputs/viz/driver", d,"_mag.png"), arclistlist, colorlist, thicknesslist, 2400, 1800)
	end

	for d in [9], f in 1:numfragments[driverHomeLocs[d], drivershift[d]]
		if value(z_smp[d,f]) > 1e-4
			println("-------- FRAG $f --------")
			println("flow = ", value(z_smp[d,f]))
			fragDesc(driverHomeLocs[d], drivershift[d],f)
		end
	end

	for i in [33], a in magarcs.A[i]
		if value(x_smp[i,a]) > 1e-4
			println("Flow = ", value(x_smp[i,a]))
			arcDesc(a)
		end
	end


	if formulation == "heterogeneous"
		println("Driver util = ", sum(sum(fragworkinghours[driverHomeLocs[d],drivershift[d],f] * value(z_magip[d,f]) for f in 1:numfragments[driverHomeLocs[d],drivershift[d]]) for d in drivers) / sum(maxweeklydriverhours for d in drivers))
	end
	=#

elseif (solutionmethod == "cg") & (formulation == "homogeneous")

	dummypath = 1
	cg_obj, rmp, x_rmp, y_rmp, z_rmp, w_rmp, cgpaths, delta, rmptime, pptime, pptime_par, totalcgpaths, cg_iter, fullcg_time = columngeneration!(orderarcs, hasdriverarcs, nocuts)
	cgip_obj, x_cgip, z_cgip, cgip_time, cgip_bound = solvejourneymodel_paths(0, opt_gap, cgpaths, delta, numeffshifts)

	timeslist1 = (mp=rmptime, pp=pptime, pppar=pptime_par, ip=0, cut=0, full=fullcg_time-pptime+pptime_par)
	writeresultsforrun(resultsfilename, 0, cg_iter, cg_obj, timeslist1, totalcgpaths, x_cgip, z_cgip)
	timeslist2 = (mp=0, pp=0, pppar=0, ip=cgip_time, cut=0, full=cgip_time)
	writeresultsforrun(resultsfilename, 1, "IP", cgip_obj, timeslist2, totalcgpaths, x_cgip, z_cgip)

elseif solutionmethod == "cg"

	dummypath = 1

	startercuts = (vars=Dict(), rhs=Dict(), coeff=Dict())
	
	cg_obj, rmp, x_rmp, y_rmp, z_rmp, w_rmp, cgpaths, delta, rmptime, pptime, pptime_par, totalcgpaths, cg_iter, knapsackcuts, cgcuttime = columngeneration!(orderarcs, hasdriverarcs, startercuts)
	cgip_obj, x_cgip, z_cgip, cgip_time, cgip_bound = solvejourneymodel_paths(0, opt_gap, cgpaths, delta, numeffshifts, knapsackcuts)

	timeslist1 = (mp=rmptime, pp=pptime, pppar=pptime_par, ip=0, cut=cgcuttime, full=0)
	writeresultsforrun(resultsfilename, 0, cg_iter, cg_obj, timeslist1, totalcgpaths, x_cgip, z_cgip)
	timeslist2 = (mp=0, pp=0, pppar=0, ip=cgip_time, cut=0, full=0)
	writeresultsforrun(resultsfilename, 1, "IP", cgip_obj, timeslist2, totalcgpaths, x_cgip, z_cgip)

end

#-----------------------------------------------------------------------------#

println("Done!")
