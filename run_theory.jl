
using JuMP, Gurobi, Random, CSV, DataFrames, Statistics, Dates, SparseArrays 

include("scripts/theoryhelper.jl")
include("scripts/visualizations/theorynetwork.jl")

data = CSV.read("data/theory.csv", DataFrame)

W = [1,2,3]
E = [4,5,6]
T = 8

demand = zeros(6, 6, T)
for row in 1:size(data)[1]
    demand[data[row,2],data[row,3],data[row,4]] = data[row,5]
end

function Tmod(a)
    return mod(a-1,T)+1
end

coordinates = [0 1 ; 0 0.5 ; 0 0 ; 3 1 ; 3 0.5 ; 3 0 ; 1 0.5; 2 0.5]

#--------------------------------------------------------------#

#Point-to-point
C = 5

#Opt
journeyscovering, journeydist = Dict(), Dict()
for i in W, j in E, t in 1:T
    journeyscovering[i,j,t] = []
    journeyscovering[j,i,t] = []
end
jindex = 1
journeyarclookup = Dict()
for i in W, j in E, t in 1:T, i2 in W, j2 in E
    push!(journeyscovering[i,j,t], jindex)
    push!(journeyscovering[j2,i2,mod(t+C-1,T)+1], jindex)
    journeydist[jindex] = sqrt(3^2 + (0.5*abs(mod(i,3)-mod(j,3)))^2) + sqrt(3^2 + (0.5*abs(mod(i2,3)-mod(j2,3)))^2) + 0.5*abs(i-i2) + 0.5*abs(j-j2)
    
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
P = [7,8]
pitstops = 1:length(union(W,E,P))
corridors = [(1,7),(2,7),(3,7),(7,8),(4,8),(5,8),(6,8)]

flow = zeros(length(union(W,E,P)), length(union(W,E,P)), T)
for i in W, j in E, t in 1:T
    flow[i,7,t] += demand[i,j,t]
    flow[7,8,Tmod(t+1)] += demand[i,j,t]
    flow[8,j,Tmod(t+2)] += demand[i,j,t] 
end
for i in E, j in W, t in 1:T
    flow[i,8,t] += demand[i,j,t]
    flow[8,7,Tmod(t+1)] += demand[i,j,t]
    flow[7,j,Tmod(t+2)] += demand[i,j,t] 
end
tripdistance = zeros(length(pitstops), length(pitstops))
for i in 1:8, j in 1:8
    tripdistance[i,j] = sqrt((coordinates[i,1] - coordinates[j,1])^2 + (coordinates[i,2] - coordinates[j,2])^2)
    tripdistance[j,i] = sqrt((coordinates[i,1] - coordinates[j,1])^2 + (coordinates[i,2] - coordinates[j,2])^2)
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
@variable(model, demandheld[i in 1:6,j in receivinglocs[i],1:T,1:T] >= 0, Int)
@objective(model, Min, sum(journeydist[j] * y[j] for j in journeys))

#Delay
@constraint(model, [i in 1:6,j in receivinglocs[i], t in 1:T], sum(demandheld[i,j,t,t2] for t2 in 1:T) == demand[i,j,t])
@constraint(model, totaldelay, sum(sum(sum(sum(periodsdelay(t,t2) * demandheld[i,j,t,t2] for t2 in 1:T) for t in 1:T) for j in receivinglocs[i]) for i in 1:6) <= maxdelay)
#@constraint(model, [i in 1:6,j in receivinglocs[i], t in 1:T], demandheld[i,j,t,t] == demand[i,j,t])
#@constraint(model, [i in 1:6,j in receivinglocs[i], t in 1:T, t2 in setdiff(1:T,t)], demandheld[i,j,t,t2] == 0)

#Flow
@constraint(model, [i in W, t in 1:T], sum(y[jrny] for jrny in journeyscovering[i,7,t]) >= sum(sum(demandheld[i,j,t2,t] for j in E) for t2 in 1:T) )
@constraint(model, [t in 1:T], sum(y[jrny] for jrny in journeyscovering[7,8,t]) >= sum(sum(sum(demandheld[i,j,Tmod(t2-1),Tmod(t-1)] for j in E) for i in W) for t2 in 1:T) )
@constraint(model, [j in E, t in 1:T], sum(y[jrny] for jrny in journeyscovering[8,j,t]) >= sum(sum(demandheld[i,j,Tmod(t2-2),Tmod(t-2)] for i in W) for t2 in 1:T) )

@constraint(model, [i in E, t in 1:T], sum(y[jrny] for jrny in journeyscovering[i,8,t]) >= sum(sum(demandheld[i,j,t2,t] for j in W) for t2 in 1:T) )
@constraint(model, [t in 1:T], sum(y[jrny] for jrny in journeyscovering[8,7,t]) >= sum(sum(sum(demandheld[i,j,Tmod(t2-1),Tmod(t-1)] for j in W) for i in E) for t2 in 1:T) )
@constraint(model, [j in W, t in 1:T], sum(y[jrny] for jrny in journeyscovering[7,j,t]) >= sum(sum(demandheld[i,j,Tmod(t2-2),Tmod(t-2)] for i in E) for t2 in 1:T) )

optimize!(model)

#println("Relay bound = ", relay_bound)
relay_delay_obj = objective_value(model)

#--------------------------------------------------------------#

#Relay bound
relay_bound = 0
for i in W, j in [7]
    r = [sum(demand[i,j,t] for j in E) for t in 1:T]
    l = [sum(demand[j,i,Tmod(t-2)] for j in E) for t in 1:T]
    minmiles_corr = 2 * tripdistance[i,j] * findminmiles(r, l)
    global relay_bound += minmiles_corr
end 
for i in [7], j in [8]
    r = [sum(sum(demand[i,j,Tmod(t-1)] for j in E) for i in W) for t in 1:T]
    l = [sum(sum(demand[i,j,Tmod(t-1)] for j in W) for i in E) for t in 1:T]
    minmiles_corr = 2 * tripdistance[i,j] * findminmiles(r, l)
    global relay_bound += minmiles_corr
end 
for i in E, j in [8]
    r = [sum(demand[i,j,t] for j in W) for t in 1:T]
    l = [sum(demand[j,i,Tmod(t-2)] for j in W) for t in 1:T]
    minmiles_corr = 2 * tripdistance[i,j] * findminmiles(r, l)
    global relay_bound += minmiles_corr
end 

#--------------------------------------------------------------#

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

if 1==1
    println("--------------------------------")
    println("Point-to-point miles = ", ptp_obj)
    #println("   Empty (ptp) = ", ptp_empties)
    println("Point-to-point bound = ", ptp_bound) 
    println("Relay bound = ", relay_bound) 
    println("Relay miles = ", relay_obj)
    #println("   Empty (relay) = ", relay_empties)
    println("Relay miles w/ delay = ", relay_delay_obj)
    #println("   Strategic delay = ", sum(sum(sum(sum(periodsdelay(t,t2) * value(demandheld[i,j,t,t2]) for t2 in 1:T) for t in 1:T) for j in receivinglocs[i]) for i in 1:6))
end