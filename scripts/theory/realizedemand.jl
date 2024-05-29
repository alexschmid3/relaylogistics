
function realizedemands(d_bar, stdev)

    numlocs = length(demandlocs)

    demand = zeros(numlocs, numlocs, T)
    for (i,j) in allpairs, t in 1:T
        #if d_bar[i,j] > 1e-4
        if demanddist == "normal"
            demand[i,j,t] = min(d_ub, max(d_lb, round(randn() * stdev_base * d_bar + d_bar, digits=0)))
        elseif demanddist == "uniform"
            demand[i,j,t] = (d_ub - d_lb) * rand() + d_lb
        end
        #end
    end

    return demand

end
