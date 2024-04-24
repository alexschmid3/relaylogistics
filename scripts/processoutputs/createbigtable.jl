
using CSV, DataFrames, Statistics

#include("createbigtable_helper.jl")
include("scripts/processoutputs/createbigtable_helper.jl")

datafile = "outputs/table2/deadlines/ex_combined_delete.csv"

fullsizes, fulltimes, fulllambdas, fullmethods = [2,3,4], [(3, 6), (4, 6), (5, 6)], [500,1000], ["Direct implementation", "Direct on LO basis", "Path-based column generation", "Single-arc generation", "Multi-arc generation", "Multi-arc generation [ENHANCED]"]

methodmap_stg = Dict("ip" => "a_ip",
                 "basisip" => "b_basisip", 
                 "cg" => "c_cg",
                 "sag" => "d_sag",
                 "mag" => "e_mag",
                 "mage" => "f_mage"
                 ) 
methodmap = Dict("a_ip" => "Direct implementation",
                 "b_basisip" => "Direct on LO basis", 
                 "c_cg" => "Path-based column generation",
                 "d_sag" => "Single-arc generation",
                 "e_mag" => "Multi-arc generation",
                 "f_mage" => "Multi-arc generation [ENHANCED]"
                 ) 

#Primary key for each instance --> "how can I identify which IP optimal solution this row corresponds to?"
#instancekey = [:instance, :lambda_delay, :lambda_drvrhrs, :horizon, :tstep, :week, :numlocs, :numorders, :numdrivers]
#Primary key for each experiment_id
#methodkey = [:experiment_id, :instance, :lambda_delay, :lambda_drvrhrs, :horizon, :tstep, :week, :numlocs, :numorders, :numdrivers, :method]
#What do you want the final table grouped by?
#tablekey = [:lambda_delay, :horizon, :tstep, :method]

#Primary key for each instance --> "how can I identify which IP optimal solution this row corresponds to?"
instancekey = [:instance, :lambda_delay, :lambda_drvrhrs, :horizon, :tstep, :week, :numlocs, :numorders, :numdrivers]
instancekey_cuts = [:instance, :lambda_delay, :lambda_drvrhrs, :horizon, :tstep, :week, :numlocs, :numorders, :numdrivers, :cuttype]
#Primary key for each experiment_id --> "over which rows should computational times be summed?"
methodkey = [:experiment_id, :instance, :lambda_delay, :lambda_drvrhrs, :horizon, :tstep, :week, :numlocs, :numorders, :numdrivers, :method]
#What do you want the final table grouped by?
tablekey = [:instance, :lambda_delay, :numdrivers, :horizon, :tstep, :method, :varsettingiterations, :variablefixingthreshold, :strongreducedcosts, :deletioniterationpercent, :deletionthreshold, :cuttype]

#----------------------------------------------------------------------------#

#Read data
df = CSV.read(datafile, DataFrame)
#df = filter(row -> row.week == 1, df)
#Get rid of missing columns 
df = replacemissingwithzeros(df, "lambda_drvrhrs")
df = filter(row -> row.lambda_delay == 500, df)

#Find the IP optimal solution for all instances
IPoptima = filter(row -> findIPopt(row.method, row.iteration) == true, df)
IPoptima = IPoptima[:, ["instance","lambda_delay","lambda_drvrhrs","horizon","tstep","week","numlocs","numorders","numdrivers","objective"]]
rename!(IPoptima,:objective => :ipopt)

#Find the LP optimal solution for all instances
#LPoptima = filter(row -> findLPopt(row.method, row.iteration) == true, df)
#LPoptima = LPoptima[:, ["instance","lambda_delay","lambda_drvrhrs","horizon","tstep","week","numlocs","numorders","numdrivers","cuttype","objective"]]
#rename!(LPoptima,:objective => :lpopt)

#LPoptima = filter(row -> row.experiment_id >= 1449, df)
LPoptima = filter(row -> row.method == "mag", df)
LPoptima = filter(row -> row.iteration != "IP", LPoptima)
LPoptima = LPoptima[:, ["instance","lambda_delay","lambda_drvrhrs","horizon","tstep","week","numlocs","numorders","numdrivers","cuttype","objective"]]
rename!(LPoptima,:objective => :lpopt)
unique!(LPoptima)

#Join the IP and LP optima onto the original table and calculate the IP and LP gaps
df = outerjoin(df, LPoptima, on = instancekey_cuts)
df = outerjoin(df, IPoptima, on = instancekey)
lpgaps = [(df[i,"objective"]-df[i,"lpopt"])/df[i,"lpopt"] for i in 1:size(df)[1]]
ipgaps = [(df[i,"objective"]-df[i,"ipopt"])/df[i,"ipopt"] for i in 1:size(df)[1]]
df[!, "lpgap"] = lpgaps
df[!, "ipgap"] = ipgaps

#Sum multiple row times per method
df_grp = groupby(df, methodkey)
df_times = combine(df_grp, [:smptime, :pptime_par, :iptime, :cuttime] => ((mp,pp,ip,cut) -> (totaltime=sum(mp)+sum(pp)+sum(ip),mptime=sum(mp), pptime=sum(pp), iptime=sum(ip), cuttime=sum(cut))) => AsTable) 

#Add the total times to the full table
df_sol = filter(row -> row.iteration == "IP", df)
df_sol = select!(df_sol, Not([:smptime, :pptime, :pptime_par, :iptime, :cuttime]));
df_table = outerjoin(df_sol, df_times, on = methodkey)

#Average over 8 weeks for each instance size
df_grp = groupby(df_table, tablekey)
df_summ = combine(df_grp, nrow, [:totaltime, :mptime, :pptime, :iptime, :cuttime, :objective, :lpgap, :ipgap] => ((tt,mp,pp,ip,cut,obj,lpgap,ipgap) -> (totaltime=mean(skipmissing(tt)),mptime=mean(skipmissing(mp)), pptime=mean(skipmissing(pp)), iptime=mean(skipmissing(ip)), cuttime=mean(skipmissing(cut)), vsopt=mean(skipmissing(ipgap)), iogap=mean(skipmissing(lpgap)))) => AsTable) 

#Sort for table
df_summ[!,"method"] = [methodmap_stg[i] for i in df_summ[!,"method"]]
#df_final = sort!(df_summ, [:horizon, order(:tstep, rev=true), :lambda_delay, :method])
df_final = sort!(df_summ, [:instance, :numdrivers, :horizon, order(:tstep, rev=true), :lambda_delay, :method, :varsettingiterations, :variablefixingthreshold, :strongreducedcosts, :deletioniterationpercent, :deletionthreshold])
#CSV.write("myfile_updated2.csv", df_final)

#Format as LaTeX table
#printbigtable(df_final, fulltimes, fulllambdas, fullmethods)
df_final[!,"method"] = [methodmap[i] for i in df_final[!,"method"]]
CSV.write("outputs/table2_deadlines_v7.csv", df_final)

