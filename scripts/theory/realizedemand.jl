
function realizedemands(d_bar, stdev)

    numlocs = length(demandlocs)

    demand = zeros(numlocs, numlocs, T)
    for (i,j) in allpairs, t in 1:T
        if d_bar[i,j] > 1e-4
            demand[i,j,t] = max(0.0, round(randn() * stdev[i,j] + d_bar[i,j], digits=0))
        end
    end

    return demand

end