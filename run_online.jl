
using JuMP, Gurobi, Random, CSV, DataFrames, Statistics, Dates, SparseArrays 

#-----------------------------------LOAD OTHER FILES------------------------------------#

include("scripts/instancegeneration/readrivigodata.jl")
include("scripts/instancegeneration/shortestpath.jl")
include("scripts/instancegeneration/constraintmatrix.jl")
include("scripts/journeybasedmodel/initializejourneymodel.jl")
include("scripts/metrics/writeresultsforrun.jl")

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
experiment_id = ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
paramsfilename = "data/table3.csv"
expparms = CSV.read(paramsfilename, DataFrame)
formulation = expparms[experiment_id, 15]  # Drivers = homogeneous, heterogeneous
ex = expparms[experiment_id, 2]		
weekstart = expparms[experiment_id, 4]
horizon = expparms[experiment_id, 5] * 24
tstep = expparms[experiment_id, 6]
driverfactor = expparms[experiment_id, 7]
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
if formulation == "heterogeneous"
	include("scripts/journeybasedmodel/solvejourneymodel.jl")
	include("scripts/journeybasedmodel/solvejourneymodel_paths.jl")
	include("scripts/multiarcgeneration/multiarcgeneration_heterogeneous.jl")
	include("scripts/columngeneration/columngeneration.jl")
	include("scripts/arcbasedmodel/solvearcbasedmodel_heterogeneous.jl")
elseif formulation == "homogeneous"
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

#Import miscellaneous functions
if maketimespacevizfiles + makespatialvizfiles + makeadvancedvizfiles + vizflag >= 1
	include("scripts/visualizations/timespacenetwork.jl")
end
include("scripts/directoryinitialization.jl")
include("scripts/helper.jl")
