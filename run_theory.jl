
using JuMP, Gurobi, Random, CSV, DataFrames, Statistics, Dates, SparseArrays 

include("scripts/theory/networkinfo.jl")
include("scripts/theory/theoryhelper.jl")
include("scripts/theory/generatedemand.jl")
include("scripts/theory/realizedemand.jl")
include("scripts/visualizations/theorynetwork.jl")

#Read experiment parameters from file
experiment_id += 1 #ifelse(length(ARGS) > 0, parse(Int, ARGS[1]), 1)
params = CSV.read("data/heatmap1.csv", DataFrame)

w = 2
h = 1/2
n = params[experiment_id, 6]
T = 1
C = 5
m = 1
K = 1

#totalflow = 500 #params[experiment_id, 5] 
stdev_base = params[experiment_id, 7] 
aggbalance = params[experiment_id, 4] 
disaggbalance = params[experiment_id, 3] 
coastbalance = params[experiment_id, 2] 
randomseedval = params[experiment_id, 8] 
demanddist = params[experiment_id, 9] 
d_lb = params[experiment_id, 10] 
d_ub = params[experiment_id, 11] 
Random.seed!(randomseedval)

W = [i for i in 1:n]
E = [i for i in n+1:2*n]

demandlocs = union(W,E)
allpairs = []
for i in W, j in E
    push!(allpairs, (i,j))
end
for i in E, j in W
    push!(allpairs, (i,j))
end

coordinates, P = networkinfo()
p1,p2 = P[1], P[2]
alllocs = union(demandlocs, P)

function Tmod(a)
    return mod(a-1,T)+1
end

#--------------------------------------------------------------#

#d_bar, stdev, actualAB, actualDB, actualCB = generatedemand(totalflow, aggbalance, disaggbalance, coastbalance) 
d_bar = (d_ub-d_lb)*0.5 + d_lb
demand = realizedemands(d_bar, stdev_base)
#demand = zeros(6,6,T)
#for t in 1:T
#    demand[1,4,t] = 1
#    demand[2,5,t] = 1
#    demand[3,6,t] = 1
#    demand[4,2,t] = 1
#    demand[5,3,t] = 1
#    demand[6,1,t] = 1
#end

#--------------------------------------------------------------#

tripdistance = zeros(length(alllocs), length(alllocs))
for i in alllocs, j in alllocs
    tripdistance[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2 + (coordinates[i,2] - coordinates[j,2])^2)
    tripdistance[j,i] = sqrt((coordinates[i,1] - coordinates[j,1])^2 + (coordinates[i,2] - coordinates[j,2])^2)
end

#--------------------------------------------------------------#

#Point-to-point
journeyscovering, journeydist = Dict(), Dict()
for i in W, j in E, t in 1:T
    journeyscovering[i,j,t] = []
    journeyscovering[j,i,t] = []
end
jindex = 1
journeyarclookup = Dict()
for i in W, j in E, t in 1:T, i2 in W, j2 in E #i2 in [i], j2 in [j]
    push!(journeyscovering[i,j,t], jindex)
    push!(journeyscovering[j2,i2,mod(t+C-1,T)+1], jindex)
    journeydist[jindex] = tripdistance[i,j] + tripdistance[i2,j2] + tripdistance[i,i2] + tripdistance[j,j2] 
    
    journeyarcs = []
    for (orig,dest,t_orig,t_dest) in [(i,j,t,Tmod(t+C)), (j,j2,Tmod(t+C),Tmod(t+C)), (j2,i2,Tmod(t+C),Tmod(t+2C)), (i2,i,Tmod(t+2*C),Tmod(t+2*C))]
        if orig != dest
            push!(journeyarcs, (orig,dest,t_orig,t_dest))
        end
    end
    journeyarclookup[jindex] = journeyarcs
    
    global jindex += 1
end 
journeys = 1:jindex-1
    
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "OutputFlag", 0)
@variable(model, y[j in journeys] >= 0, Int)
@objective(model, Min, sum(journeydist[j] * y[j] for j in journeys))
@constraint(model, [i in W, j in E, t in 1:T], sum(y[jrny] for jrny in journeyscovering[i,j,t]) >= demand[i,j,t])
@constraint(model, [i in E, j in W, t in 1:T], sum(y[jrny] for jrny in journeyscovering[i,j,t]) >= demand[i,j,t])

optimize!(model)

ptp_obj = objective_value(model)

ptp_empties = 0 
for i in W, j in E, t in 1:T
    global ptp_empties += (sum(value(y[jrny]) for jrny in journeyscovering[i,j,t]) - demand[i,j,t]) * sqrt((coordinates[i,1] - coordinates[j,1])^2 + (coordinates[i,2] - coordinates[j,2])^2)
end
for i in E, j in W, t in 1:T
    global ptp_empties += (sum(value(y[jrny]) for jrny in journeyscovering[i,j,t]) - demand[i,j,t]) * sqrt((coordinates[i,1] - coordinates[j,1])^2 + (coordinates[i,2] - coordinates[j,2])^2)
end

theorynetwork("pointtopoint.png", y, 2000, 2000)

#--------------------------------------------------------------#

#Relay
pitstops = 1:length(union(W,E,P))
corridors = [] 
for i in W
    push!(corridors, (i,P[1]))
end
push!(corridors, (P[1], P[2]))
for i in E
    push!(corridors, (i,P[2]))
end

flow = zeros(length(union(W,E,P)), length(union(W,E,P)), T)
for i in W, j in E, t in 1:T
    flow[i,p1,t] += demand[i,j,t]
    flow[p1,p2,Tmod(t+1)] += demand[i,j,t]
    flow[p2,j,Tmod(t+2)] += demand[i,j,t] 
end
for i in E, j in W, t in 1:T
    flow[i,p2,t] += demand[i,j,t]
    flow[p2,p1,Tmod(t+1)] += demand[i,j,t]
    flow[p1,j,Tmod(t+2)] += demand[i,j,t] 
end

journeyscovering, journeydist = Dict(), Dict()
for (i,j) in corridors, t in 1:T
    journeyscovering[i,j,t] = []
    journeyscovering[j,i,t] = []
end
jindex = 1
for (i,j) in corridors, t in 1:T
    push!(journeyscovering[i,j,t], jindex)
    push!(journeyscovering[j,i,Tmod(t+1)], jindex)
    journeydist[jindex] = tripdistance[i,j] + tripdistance[j,i]
    global jindex += 1
end 
for (j,i) in corridors, t in 1:T
    push!(journeyscovering[i,j,t], jindex)
    push!(journeyscovering[j,i,Tmod(t+1)], jindex)
    journeydist[jindex] = tripdistance[i,j] + tripdistance[j,i]
    global jindex += 1
end 
journeys = 1:jindex-1
    
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "OutputFlag", 0)
@variable(model, y[j in journeys] >= 0, Int)
@objective(model, Min, sum(journeydist[j] * y[j] for j in journeys))
@constraint(model, [(i,j) in corridors, t in 1:T], sum(y[jrny] for jrny in journeyscovering[i,j,t]) >= flow[i,j,t])
@constraint(model, [(j,i) in corridors, t in 1:T], sum(y[jrny] for jrny in journeyscovering[i,j,t]) >= flow[i,j,t])

optimize!(model)

#println("Relay bound = ", relay_bound)
relay_obj = objective_value(model)

relay_empties = 0
for (i,j) in corridors, t in 1:T
    global relay_empties += (sum(value(y[jrny]) for jrny in journeyscovering[i,j,t]) - flow[i,j,t]) * tripdistance[i,j]
end
for (j,i) in corridors, t in 1:T
    global relay_empties += (sum(value(y[jrny]) for jrny in journeyscovering[i,j,t]) - flow[i,j,t]) * tripdistance[i,j]
end

#--------------------------------------------------------------#
#=
function periodsdelay(t1,t2)
    if t2 >= t1
        return t2 - t1
    else
        return t2 + T - t1
    end
end

#Relay w/ delay
maxdelay = 2 * sum(demand)
receivinglocs = Dict()
for i in W
    receivinglocs[i] = E
end
for i in E
    receivinglocs[i] = W
end
    
model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "OutputFlag", 0)
@variable(model, y[j in journeys] >= 0, Int)
@variable(model, demandheld[i in demandlocs, j in receivinglocs[i],1:T,1:T] >= 0, Int)
@objective(model, Min, sum(journeydist[j] * y[j] for j in journeys))

#Delay
@constraint(model, [i in demandlocs, j in receivinglocs[i], t in 1:T], sum(demandheld[i,j,t,t2] for t2 in 1:T) == demand[i,j,t])
@constraint(model, totaldelay, sum(sum(sum(sum(periodsdelay(t,t2) * demandheld[i,j,t,t2] for t2 in 1:T) for t in 1:T) for j in receivinglocs[i]) for i in demandlocs) <= maxdelay)
#@constraint(model, [i in demandlocs,j in receivinglocs[i], t in 1:T], demandheld[i,j,t,t] == demand[i,j,t])
#@constraint(model, [i in demandlocs,j in receivinglocs[i], t in 1:T, t2 in setdiff(1:T,t)], demandheld[i,j,t,t2] == 0)

#Flow
@constraint(model, [i in W, t in 1:T], sum(y[jrny] for jrny in journeyscovering[i,p1,t]) >= sum(sum(demandheld[i,j,t2,t] for j in E) for t2 in 1:T) )
@constraint(model, [t in 1:T], sum(y[jrny] for jrny in journeyscovering[p1,p2,t]) >= sum(sum(sum(demandheld[i,j,Tmod(t2-1),Tmod(t-1)] for j in E) for i in W) for t2 in 1:T) )
@constraint(model, [j in E, t in 1:T], sum(y[jrny] for jrny in journeyscovering[p2,j,t]) >= sum(sum(demandheld[i,j,Tmod(t2-2),Tmod(t-2)] for i in W) for t2 in 1:T) )

@constraint(model, [i in E, t in 1:T], sum(y[jrny] for jrny in journeyscovering[i,p2,t]) >= sum(sum(demandheld[i,j,t2,t] for j in W) for t2 in 1:T) )
@constraint(model, [t in 1:T], sum(y[jrny] for jrny in journeyscovering[p2,p1,t]) >= sum(sum(sum(demandheld[i,j,Tmod(t2-1),Tmod(t-1)] for j in W) for i in E) for t2 in 1:T) )
@constraint(model, [j in W, t in 1:T], sum(y[jrny] for jrny in journeyscovering[p1,j,t]) >= sum(sum(demandheld[i,j,Tmod(t2-2),Tmod(t-2)] for i in E) for t2 in 1:T) )

optimize!(model)
relay_delay_obj = objective_value(model)
=#
#--------------------------------------------------------------#
#=
#Relay bound
relay_bound = 0
for i in W, j in [p1]
    r = [sum(demand[i,j,t] for j in E) for t in 1:T]
    l = [sum(demand[j,i,Tmod(t-2)] for j in E) for t in 1:T]
    minmiles_corr = 2 * tripdistance[i,j] * findminmiles(r, l)
    global relay_bound += minmiles_corr
end 
for i in [p1], j in [p2]
    r = [sum(sum(demand[i,j,Tmod(t-1)] for j in E) for i in W) for t in 1:T]
    l = [sum(sum(demand[i,j,Tmod(t-1)] for j in W) for i in E) for t in 1:T]
    minmiles_corr = 2 * tripdistance[i,j] * findminmiles(r, l)
    global relay_bound += minmiles_corr
end 
for i in E, j in [p2]
    r = [sum(demand[i,j,t] for j in W) for t in 1:T]
    l = [sum(demand[j,i,Tmod(t-2)] for j in W) for t in 1:T]
    minmiles_corr = 2 * tripdistance[i,j] * findminmiles(r, l)
    global relay_bound += minmiles_corr
end 
=#
#--------------------------------------------------------------#
#=
#Ptp bound 
we_bound, ew_bound = 0, 0
for i in W, j in E, t in 1:T
    totalneed = demand[i,j,t]
    options = []
    for i2 in W, j2 in E, t2 in [Tmod(t+C), Tmod(t-C)]
        totaldistance = tripdistance[j,j2] + tripdistance[i2,i]
        maxmatch = demand[j2,i2,t2]
        push!(options, ((j2,i2,t2), totaldistance, maxmatch))
    end

    sortedoptions = sort(options, by=x->x[2])
    push!(sortedoptions, ("empty",  tripdistance[i,j], totalneed))

    filledneed, totalmiles = 0, 0
    for opt in sortedoptions
        need = min(totalneed - filledneed, opt[3])
        filledneed += need
        totalmiles += need * opt[2]
    end

    global we_bound += totalmiles
end
for i in E, j in W, t in 1:T
    totalneed = demand[i,j,t]
    options = []
    for i2 in E, j2 in W, t2 in [Tmod(t+C), Tmod(t-C)]
        totaldistance = tripdistance[j,j2] + tripdistance[i2,i]
        maxmatch = demand[j2,i2,t2]
        push!(options, ((j2,i2,t2), totaldistance, maxmatch))
    end

    sortedoptions = sort(options, by=x->x[2])
    push!(sortedoptions, ("empty",  tripdistance[i,j], totalneed))

    filledneed, totalmiles = 0, 0
    for opt in sortedoptions
        need = min(totalneed - filledneed, opt[3])
        filledneed += need
        totalmiles += need * opt[2]
    end

    global ew_bound += totalmiles
end
ptp_bound = max(we_bound, ew_bound) + sum(sum(sum(tripdistance[i,j]*demand[i,j,t] for t in 1:T) for j in E) for i in W) + sum(sum(sum(tripdistance[i,j]*demand[i,j,t] for t in 1:T) for j in W) for i in E) 

#--------------------------------------------------------------#

#PTP bound 2
α = 0.5
M = 2*sum(sum(sum(demand[i,j,t] * tripdistance[i,j] for t in 1:T) for j in demandlocs) for i in demandlocs)

#Bipartite graph - vertices
vertexid, vertexdesc = Dict(), Dict()
U, V = [], []
vindex = 1
for i in W, j in E, t in 1:T, k in 1:demand[i,j,t]
    vertexdesc[vindex] = (i,j,t,k)
    vertexid[i,j,t,k] = vindex
    push!(U, vindex)
    global vindex += 1
end
for i in E, j in W, t in 1:T, k in 1:demand[i,j,t]
    vertexdesc[vindex] = (i,j,t,k)
    vertexid[i,j,t,k] = vindex
    push!(V, vindex)
    global vindex += 1
end
vertices = 1:vindex-1

#Bipartite graph - edges
edges = []
milessaved, pairsof = Dict(), Dict()
for i in vertices
    pairsof[i] = []
end
for i in W, j in E, t in 1:T, k in 1:demand[i,j,t], i2 in W, j2 in E, k2 in 1:demand[j2,i2,Tmod(t+C)]
    vert1 = vertexid[i,j,t,k]
    vert2 = vertexid[j2,i2,Tmod(t+C),k2]
    push!(edges, (vert1, vert2))
    push!(pairsof[vert1], vert2)
    push!(pairsof[vert2], vert1)
    milessaved[vert1, vert2] = tripdistance[i,j] + tripdistance[j2,i2] - tripdistance[i,i2] - tripdistance[j,j2]
end
for i in W, j in E, t in 1:T, k in 1:demand[i,j,t], i2 in W, j2 in E, k2 in 1:demand[j2,i2,Tmod(t-C)]
    vert1 = vertexid[i,j,t,k]
    vert2 = vertexid[j2,i2,Tmod(t-C),k2]
    push!(edges, (vert1, vert2))
    push!(pairsof[vert1], vert2)
    push!(pairsof[vert2], vert1)
    milessaved[vert1, vert2] = tripdistance[i,j] + tripdistance[j2,i2] - tripdistance[i,i2] - tripdistance[j,j2]
end

#Calculation of bound
ptp_bound2 = 0
for α in 0:0.1:1
    pvar, qvar = Dict(), Dict()
    for u in U
        pvar[u] = maximum([α * milessaved[u,v] for v in pairsof[u]])
    end
    for v in V
        qvar[v] = maximum([(1-α) * milessaved[u,v] for u in pairsof[v]])
    end
    global ptp_bound2 = max(ptp_bound2, M - sum(pvar[u] for u in U) - sum(qvar[v] for v in V))
end
=#
#--------------------------------------------------------------#
#=
using Hungarian

weights = Array{Union{Missing, Float64}}(undef,length(U), length(V)) 
shift = length(U)
for (u,vind) in edges
    v = vind-shift 
    weights[u,v] = milessaved[u,vind]
end    
assignments, cost = hungarian(weights)

ptp_exact = M - cost
=#
#--------------------------------------------------------------#
#=
#PTP simple (T=C)
ptp_bound3 = sum(sum(sum(tripdistance[i,j]*demand[i,j,t] for t in 1:T) for j in E) for i in W) + sum(sum(sum(tripdistance[i,j]*demand[i,j,t] for t in 1:T) for j in W) for i in E) 

disagg = sum(sum(sum(abs(demand[i,j,t] - demand[j,i,t]) for i in W) for j in E) for t in 1:T)
coastal = sum(abs(sum(sum(demand[i,j,t] for j in E) for i in W) - sum(sum(demand[i,j,t] for j in W) for i in E)) for t in 1:T)

ptp_bound3 += (3*w - (h)/(n-1)) * coastal + (h)/(n-1) * disagg

#Relay simple
detourdistance = Dict()
for i in W, j in E
    detourdistance[i,j] = tripdistance[i,P[1]] + tripdistance[P[1],P[2]] + tripdistance[P[2],j]
    detourdistance[j,i] = tripdistance[i,P[1]] + tripdistance[P[1],P[2]] + tripdistance[P[2],j]
end
relay_bound2 = sum(sum(sum(detourdistance[i,j]*demand[i,j,t] for t in 1:T) for j in E) for i in W) + sum(sum(sum(tripdistance[i,j]*demand[i,j,t] for t in 1:T) for j in W) for i in E) 

agg = sum(sum(abs(sum(demand[i,j,t] - demand[j,i,t] for j in E)) for i in W)  for t in 1:T) + sum(sum(abs(sum(demand[i,j,t] - demand[j,i,t] for j in W)) for i in E) for t in 1:T)

relay_bound2 += w * coastal + sqrt(w^2 + h^2) * agg
=#
#--------------------------------------------------------------#

if 1==1
    println("--------------------------------")
    #println("Maximum possible miles = ", M)
    println("Point-to-point miles = ", ptp_obj)
    #println("Point-to-point bound (exact) = ", ptp_exact)
    #println("   Empty (ptp) = ", ptp_empties)
    #println("Point-to-point bound (simple) = ", ptp_bound3) 
    #println("Point-to-point bound (new) = ", ptp_bound2) 
    #println("Point-to-point bound = ", ptp_bound) 
    #println("Relay bound = ", relay_bound) 
    #println("Relay bound (simple) = ", relay_bound2) 
    println("Relay miles = ", relay_obj)
    #println("   Empty (relay) = ", relay_empties)
    #println("Relay miles w/ delay = ", relay_delay_obj)
    #println("   Strategic delay = ", sum(sum(sum(sum(periodsdelay(t,t2) * value(demandheld[i,j,t,t2]) for t2 in 1:T) for t in 1:T) for j in receivinglocs[i]) for i in demandlocs))
    println("Minimum possible miles = ", sum(sum(sum(tripdistance[i,j]*demand[i,j,t] for t in 1:T) for j in E) for i in W) + sum(sum(sum(tripdistance[i,j]*demand[i,j,t] for t in 1:T) for j in W) for i in E) )
end

#--------------------------------------------------------------#

#df = DataFrame(ab=[actualAB], db=[actualDB], cd=[actualCB], flow=[totalflow], stdev=[stdev_base], relay=[relay_delay_obj], ptp=[ptp_obj], relay_bound=[relay_bound2], ptp_bound=[ptp_bound3])
df = DataFrame(experiment_id=[experiment_id], seed=[randomseedval], ab=[aggbalance], db=[disaggbalance], cd=[coastbalance], n=[n], stdev=[stdev_base], relay=[relay_obj], ptp=[ptp_obj], improvement=[(relay_obj-ptp_obj)/ptp_obj])
#CSV.write(string("outputs/heatmapdata/exp", experiment_id,".csv"), df)
CSV.write(string("outputs/heatmapdata/heatmap1_repos_outputs.csv"), df, append=true)
