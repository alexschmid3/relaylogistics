
function missinggap(solutionmethod, ipgap)::Bool
    interesting_name = solutionmethod == "mag"
    interesting_iteration = ismissing(ipgap)
    interesting_name && interesting_iteration
end

#----------------------------------------------------------------------------------------#

function findIPopt(solutionmethod, iteration)::Bool
    interesting_name = solutionmethod == "ip"
    interesting_iteration = "IP" == iteration
    interesting_name && interesting_iteration
end

#----------------------------------------------------------------------------------------#

function findLPopt(solutionmethod, iteration)::Bool
    interesting_name = solutionmethod == "basisip"
    interesting_iteration = "LP" == iteration
    interesting_name && interesting_iteration
end

#----------------------------------------------------------------------------------------#

function replacemissingwithzeros(df, columnname)
    
    newcolumn = []
    for i in 1:size(df)[1]
        if ismissing(df[i,columnname])
            push!(newcolumn, 0)
        else
            push!(newcolumn, df[i,columnname])
        end
    end

    df[!,"lambda_drvrhrs"] = newcolumn

    return df
end

#----------------------------------------------------------------------------------------#

function convertorhyphen(newtype, value)
    if ismissing(value)
        return "---"
    else
        return convert(newtype, value)
    end
end

#----------------------------------------------------------------------------------------#

function printbigtable(df_final, fulltimes, fulllambdas, fullmethods)

    println("Horizon & Discretization & \$ \\lambda\$ & Method & Total & MP & PP & IO  & Success & vs. OPT  & IO gap \\")
    
    currhorizon, currtstep, currlambda = 0, 0, 0

    for (hor,disc) in fulltimes, lamb in fulllambdas, meth in fullmethods
        println("$hor, $disc, $lamb, $meth")
    end    

    (hor, disc) = fulltimes[1] 
    lamb, meth = fulllambdas[1], fullmethods[1]
    for row in eachrow(df_final)
        println(convert(Int, row["horizon"]/24) == hor) 
        println(convert(Int, row["tstep"]) == disc)  
        println(convert(Int, row["lambda_delay"]) == lamb) 
        println(methodmap[row["method"]] == meth)
    end
    
    for (hor,disc) in fulltimes, lamb in fulllambdas, meth in fullmethods

        rowexists_flag = 0
        for row in eachrow(df_final)

            if (convert(Int, row["horizon"]/24) == hor) & (convert(Int, row["tstep"]) == disc) & (convert(Int, row["lambda_delay"]) == lamb) & (methodmap[row["method"]] == meth)

                if row["horizon"] != currhorizon
                    println("\\midrule")
                    print(convert(Int, row["horizon"]/24), " days")
                    print(" & ")
                    print(convert(Int, row["tstep"]), " hours")
                    print(" & ")
                    print(convert(Int, row["lambda_delay"]))
                    currhorizon = row["horizon"]
                    currtstep = row["tstep"]
                    currlambda = row["lambda_delay"]
                elseif row["tstep"] != currtstep
                    println("\\cmidrule(lr){2-11}")
                    print(" & ")
                    print(convert(Int, row["tstep"]), " hours")
                    print(" & ")
                    print(convert(Int, row["lambda_delay"]))
                    currtstep = row["tstep"]
                    currlambda = row["lambda_delay"]
                elseif row["lambda_delay"] != currlambda
                    println("\\cmidrule(lr){3-11}")
                    print(" & ")
                    print(" & ")
                    print(convert(Int, row["lambda_delay"]))
                    currlambda = row["lambda_delay"]
                else
                    print(" &  & ")
                end
                print(" & ")
                if row["method"] == "e_mag"
                    print("\\textbf{", methodmap[row["method"]], "}")
                    print(" & ")
                    print("\\textbf{", convert(Int, round(row["totaltime"]/60, digits=0)), "}")
                    print(" & ")
                    print("\\textbf{", convert(Int, round(row["mptime"]/60, digits=0)), "}")
                    print(" & ")
                    print("\\textbf{", convert(Int, round(row["pptime"]/60, digits=0)), "}")
                    print(" & ")
                    print("\\textbf{", convert(Int, round(row["iptime"]/60, digits=0)), "}")
                    print(" & ")
                    print("\\textbf{", convert(Int, row["nrow"]), "/8", "}")
                    print(" & ")
                    print("\\textbf{", convertorhyphen(Float64, round(100*row["vsopt"], digits=2)), "\\%", "}")
                    print(" & ")
                    print("\\textbf{", convertorhyphen(Float64, round(100*row["iogap"], digits=2)), "\\%", "}")
                else 
                    print(methodmap[row["method"]])
                    print(" & ")
                    print(convert(Int, round(row["totaltime"]/60, digits=0)))
                    print(" & ")
                    print(convert(Int, round(row["mptime"]/60, digits=0)))
                    print(" & ")
                    print(convert(Int, round(row["pptime"]/60, digits=0)))
                    print(" & ")
                    print(convert(Int, round(row["iptime"]/60, digits=0)))
                    print(" & ")
                    print(convert(Int, row["nrow"]), "/8")
                    print(" & ")
                    print(convertorhyphen(Float64, round(100*row["vsopt"], digits=2)), "\\%")
                    print(" & ")
                    print(convertorhyphen(Float64, round(100*row["iogap"], digits=2)), "\\%")
                end
                println(" \\\\ ")

                rowexists_flag = 1
                break

            end
        end

        if rowexists_flag == 0

            #=println("    --> ", hor*24, " vs. ", currhorizon)
            println("    --> ", disc, " vs. ", currtstep)
            println("    --> ", lamb, " vs. ", currlambda)

            println("    --> ", typeof(hor*24), " vs. ", typeof(currhorizon))
            println("    --> ", typeof(disc), " vs. ", typeof(currtstep))
            println("    --> ", typeof(lamb), " vs. ", typeof(currlambda))=#

            if hor*24 != currhorizon
                println("\\midrule")
                println(hor, " & ", disc, " & ", lamb, " & ", meth, " & --- & --- & --- & ---  & ?/8 & --- & --- \\\\")
            elseif disc != currtstep
                println("\\cmidrule(lr){2-11}")
                println("  & ", disc, " & ", lamb, " & ", meth, " & --- & --- & --- & ---  & ?/8 & --- & --- \\\\")
            elseif lamb != currlambda
                println("\\cmidrule(lr){3-11}")
                println("  &  & ", lamb, " & ", meth, " & --- & --- & --- & ---  & ?/8 & --- & --- \\\\")
            end

            println(" &  &  & ", meth, " & --- & --- & --- & ---  & ?/8 & --- & --- \\\\")
            
            currhorizon = hor*24
            currtstep = disc
            currlambda = lamb
        end

    end

end

    
