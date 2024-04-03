
using CSV

function generatedemand(totalflow, aggbalance, disaggbalance, coastbalance) 

    locations = union(W,E)

	model = Model(Gurobi.Optimizer)
	set_optimizer_attribute(model, "OutputFlag", 0)

	@variable(model, 0 <= x[i = locations, j = locations] <= round(rand()*100, digits=0) + 200)
	@variable(model, z[i = locations] >= 0)
	@variable(model, z_int[i = locations], Bin)
	@variable(model, q[i = 1:length(locations)-1, j = i+1:length(locations)] >= 0)
	@variable(model, q_int[i = 1:length(locations)-1, j = i+1:length(locations)], Bin)
	@variable(model, p >= 0)
	@variable(model, p_int, Bin)
	@variable(model, AB)
	@variable(model, DB)
	@variable(model, CB)

	@objective(model, Min, AB + DB + CB ) #- sum(sum(rand() * x[i,j] for j in locations) for i in locations))

	@constraint(model, aggflow1[i in locations], z[i] >= sum(x[i,j] for j in locations) - sum(x[j,i] for j in locations) - totalflow * z_int[i] )
	@constraint(model, aggflow2[i in locations], z[i] <= sum(x[i,j] for j in locations) - sum(x[j,i] for j in locations) + totalflow * z_int[i] )
	@constraint(model, aggflow3[i in locations], z[i] >= sum(x[j,i] for j in locations) - sum(x[i,j] for j in locations) - totalflow * (1 - z_int[i]) )
	@constraint(model, aggflow4[i in locations], z[i] <= sum(x[j,i] for j in locations) - sum(x[i,j] for j in locations) + totalflow * (1 - z_int[i]) )

	@constraint(model, diaggflow1[i = 1:length(locations)-1, j = i+1:length(locations)], q[i,j] >= x[i,j] - x[j,i] - totalflow * q_int[i,j] )
	@constraint(model, diaggflow2[i = 1:length(locations)-1, j = i+1:length(locations)], q[i,j] <= x[i,j] - x[j,i] + totalflow * q_int[i,j] )
	@constraint(model, diaggflow3[i = 1:length(locations)-1, j = i+1:length(locations)], q[i,j] >= x[j,i] - x[i,j] - totalflow * (1 - q_int[i,j]) )
	@constraint(model, diaggflow4[i = 1:length(locations)-1, j = i+1:length(locations)], q[i,j] <= x[j,i] - x[i,j] + totalflow * (1 - q_int[i,j]) )

	@constraint(model, coastflow1, p >= sum(sum(x[i,j] for j in E) for i in W) - sum(sum(x[j,i] for j in E) for i in W) - totalflow * p_int )
	@constraint(model, coastflow2, p <= sum(sum(x[i,j] for j in E) for i in W) - sum(sum(x[j,i] for j in E) for i in W) + totalflow * p_int )
	@constraint(model, coastflow3, p >= sum(sum(x[j,i] for j in E) for i in W) - sum(sum(x[i,j] for j in E) for i in W) - totalflow * (1 - p_int) )
	@constraint(model, coastflow4, p <= sum(sum(x[j,i] for j in E) for i in W) - sum(sum(x[i,j] for j in E) for i in W) + totalflow * (1 - p_int) )

	@constraint(model, westflow[i in W, j in W], x[i,j] == 0)
	@constraint(model, eastflow[i in E, j in E], x[i,j] == 0)

	@constraint(model, sum(sum(x[i,j] for j in locations) for i in locations) == totalflow)

	@constraint(model, aggbalancecalc1, AB >= sum(z[i] for i in locations) - (1-aggbalance) * totalflow)
	@constraint(model, aggbalancecalc2, AB >= (1-aggbalance) * totalflow - sum(z[i] for i in locations))
	@constraint(model, disaggbalancecalc1, DB >= sum(sum(q[i,j] for j = i+1:length(locations)) for i in 1:length(locations)-1) - (1-disaggbalance) * totalflow)
	@constraint(model, disaggbalancecalc2, DB >= (1-disaggbalance) * totalflow - sum(sum(q[i,j] for j = i+1:length(locations)) for i in 1:length(locations)-1))
	@constraint(model, coastbalancecalc1, CB >= p - (1-coastbalance) * totalflow)
	@constraint(model, coastbalancecalc2, CB >= (1-coastbalance) * totalflow - p)

	#Spread the demand
	@constraint(model, [(i,j) in allpairs], x[i,j] >= totalflow / (length(allpairs) * 2) * disaggbalance )

	#---------------------------------------------#

	optimize!(model)

	println(objective_value(model))

	#---------------------------------------------#

	demandflow = Array(value.(x))

	actualAB = round(1 - sum(abs(sum(value(x[i,j]) for j in locations) - sum(value(x[j,i]) for j in locations)) for i in locations) / totalflow, digits=3)
	actualDB = round(1 - sum(sum(abs(value(x[i,j]) - value(x[j,i])) for j in i+1:length(locations)) for i in 1:length(locations)-1) / totalflow, digits=3)
	actualCB = round(1 - abs(sum(sum(value(x[i,j]) for j in E) for i in W) - sum(sum(value(x[j,i]) for j in E) for i in W)) / totalflow, digits=3)
	println("Aggreg balance = ", actualAB)
	println("Disagg balance = ", actualDB)
	println("Coastl balance = ", actualCB)

	#CSV.write(string("data/withcoast/scenario_", actualAB, "_", actualDB, "_", actualCB, ".csv"), Tables.table(demandflow))
	#CSV.write(demandfilename, Tables.table(demandflow))

	#---------------------------------------------#

    d_bar = value.(x)
	stdev = stdev_base .* d_bar

    #---------------------------------------------#

	return d_bar, stdev, actualAB, actualDB, actualCB

end

