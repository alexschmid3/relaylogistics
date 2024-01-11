
using CSV, DataFrames
#include("createbigtable_helper.jl")
include("scripts/processoutputs/createbigtable_helper.jl")

datafile = "outputs/bigtable_orig/ex_combined.csv"

fulltimes, fulllambdas, fullmethods = [(5, 6), (5, 3), (7, 6)], [100,500,1000], ["Direct implementation", "Direct on LO basis", "Path-based column generation", "Single-arc generation", "Multi-arc generation"]

methodmap_stg = Dict("ip" => "a_ip",
                 "basisip" => "b_basisip", 
                 "cg" => "c_cg",
                 "sag" => "d_sag",
                 "mag" => "e_mag",
                 ) 
methodmap = Dict("a_ip" => "Direct implementation",
                 "b_basisip" => "Direct on LO basis", 
                 "c_cg" => "Path-based column generation",
                 "d_sag" => "Single-arc generation",
                 "e_mag" => "Multi-arc generation",
                 ) 

instancekey = [:instance, :lambda_delay, :lambda_drvrhrs, :horizon, :tstep, :week, :numlocs, :numorders, :numdrivers, :variablefixingthreshold]
methodkey = [:experiment_id, :instance, :lambda_delay, :lambda_drvrhrs, :horizon, :tstep, :week, :numlocs, :numorders, :numdrivers, :variablefixingthreshold, :method]
tablekey = [:lambda_delay, :horizon, :tstep, :method]

#Read data
df = CSV.read(datafile, DataFrame)

#Get rid of missing columns 
df = replacemissingwithzeros(df, "lambda_drvrhrs")

#Find the IP optimal solution for all instances
IPoptima = filter(row -> findIPopt(row.method, row.iteration) == true, df)
IPoptima = IPoptima[:, ["instance","lambda_delay","lambda_drvrhrs","horizon","tstep","week","numlocs","numorders","numdrivers","variablefixingthreshold","objective"]]
rename!(IPoptima,:objective => :ipopt)

#Find the LP optimal solution for all instances
LPoptima = filter(row -> findLPopt(row.method, row.iteration) == true, df)
LPoptima = LPoptima[:, ["instance","lambda_delay","lambda_drvrhrs","horizon","tstep","week","numlocs","numorders","numdrivers","variablefixingthreshold","objective"]]
rename!(LPoptima,:objective => :lpopt)

#Join the IP and LP optima onto the original table and calculate the IP and LP gaps
df = outerjoin(df, LPoptima, on = instancekey)
df = outerjoin(df, IPoptima, on = instancekey)
lpgaps = [(df[i,"objective"]-df[i,"lpopt"])/df[i,"lpopt"] for i in 1:size(df)[1]]
ipgaps = [(df[i,"objective"]-df[i,"ipopt"])/df[i,"ipopt"] for i in 1:size(df)[1]]
df[!, "lpgap"] = lpgaps
df[!, "ipgap"] = ipgaps

#Sum multiple row times per method
df_grp = groupby(df, methodkey)
df_times = combine(df_grp, [:smptime, :pptime_par, :iptime] => ((mp,pp,ip) -> (totaltime=sum(mp)+sum(pp)+sum(ip),mptime=sum(mp), pptime=sum(pp), iptime=sum(ip))) => AsTable) 

#Add the total times to the full table
df_sol = filter(row -> row.iteration == "IP", df)
df_sol = select!(df_sol, Not([:smptime, :pptime, :pptime_par, :iptime]));
df_table = outerjoin(df_sol, df_times, on = methodkey)
CSV.write("myfile.csv", df_table)

#Average over 8 weeks for each instance size
df_grp = groupby(df_table, tablekey)
df_summ = combine(df_grp, nrow, [:totaltime, :mptime, :pptime, :iptime, :objective, :lpgap, :ipgap] => ((tt,mp,pp,ip,obj,lpgap,ipgap) -> (totaltime=mean(skipmissing(tt)),mptime=mean(skipmissing(mp)), pptime=mean(skipmissing(pp)), iptime=mean(skipmissing(ip)), vsopt=mean(skipmissing(ipgap)), iogap=mean(skipmissing(lpgap)))) => AsTable) 

#Sort for table
df_summ[!,"method"] = [methodmap_stg[i] for i in df_summ[!,"method"]]
df_final = sort!(df_summ, [:horizon, order(:tstep, rev=true), :lambda_delay, :method])

#Format as LaTeX table
printbigtable(df_final, fulltimes, fulllambdas, fullmethods)
