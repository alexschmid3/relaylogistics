
using Plots, StatsPlots

T = 10
d, σ = 10, 5

function Tmod(a)
    return mod(a-1,T)+1
end

allrealizations = []

for iter in 1:10000
    demand = [max(0,randn()*σ+d) for t in 1:T]

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)
    @variable(model, x[1:T,1:T] >= 0)
    @variable(model, z >= 0)
    @objective(model, Min, z)
    @constraint(model, [t1 in 1:T], sum(x[t1,t2] for t2 in 1:T) == demand[t1])
    @constraint(model, [t1 in 1:T, t2 in setdiff(1:T, [t1,Tmod(t1+1),Tmod(t1+2)])], x[t1,t2] == 0)
    @constraint(model, [t1 in 1:T, t2 in [Tmod(t1+1),Tmod(t1+2)]], x[t1,t2] <= max(0,demand[t1]))
    @constraint(model, [t1 in 1:T, t2 in setdiff(1:T, t1)], x[t1,t2] >= 0)
    @constraint(model, [t2 in 1:T], z >= sum(x[t1,t2] for t1 in 1:T) - sum(demand))
    @constraint(model, [t2 in 1:T], z >= sum(demand) - sum(x[t1,t2] for t1 in 1:T))

    optimize!(model)
    shifted = [sum(value(x[t1,t2]) for t1 in 1:T) for t2 in 1:T]

    for item in shifted
        push!(allrealizations, item)
    end

    #=d_bar = mean(shifted)
    if maximum([item-d_bar for item in shifted]) < 1e-5
        println("good")
    else
        flag1 = maximum(demand) <= 3*d_bar 
        flag2 = true
        for t in 1:T
            if d_bar - demand[t] >= demand[Tmod(t-1)] + demand[Tmod(t-1)]
                flag2 = false
                break
            end
        end
        if (flag1) & (flag2)
            println("------------------------")
            println("Condition 1 = ", )
            
            println("Condition 2 = ", flag2)
            println("------------------------")
            println("Max = ", maximum(demand))
            println("Cap = ", 3*d_bar)
            println("Diff = ", maximum(demand) - 3*d_bar)
            println("------------------------")
            println(demand)
            println(shifted)
            println("------------------------")
        end
    end=#
end

histogram(allrealizations, bins = -5:0.5:25)

truedist = [randn()*σ + d for item in 1:1000000]
histogram(truedist, bins = -5:0.5:25)


