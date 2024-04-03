

function findmaxedgeweight_oddLs(r, l)
    maxnum, maxindex, maxside = 0, -1, "na"
    for i in 1:2:T
        if l[i] > maxnum
            maxnum = l[i]
            maxindex = i
            maxside = "L"
        end
    end
    for i in 2:2:T
        if r[i] > maxnum
            maxnum = r[i]
            maxindex = i
            maxside = "R"
        end
    end
    return maxnum, maxindex, maxside
end

function findmaxedgeweight_oddRs(r, l)
    maxnum, maxindex, maxside = 0, -1, "na"
    for i in 1:2:T
        if r[i] > maxnum
            maxnum = r[i]
            maxindex = i
            maxside = "R"
        end
    end
    for i in 2:2:T
        if l[i] > maxnum
            maxnum = l[i]
            maxindex = i
            maxside = "L"
        end
    end
    return maxnum, maxindex, maxside
end

function findminmiles(r, l)

    a, b = zeros(T), zeros(T)

    maxnum, maxindex, maxside = findmaxedgeweight_oddRs(r, l)
    a_o, b_o = zeros(T), zeros(T)
    minmiles_odd = 1e10
    if maxside == "R"
        for shift in 0:maxnum
            a_o[maxindex] = r[maxindex] - shift
            for t in [Tmod(t2) for t2 in maxindex+2:2:maxindex+T-1]
                b_o[Tmod(t-1)] = max(0, l[Tmod(t-1)]-a_o[Tmod(t-2)])
                a_o[t] = max(0, r[t]-b_o[Tmod(t-1)])
            end
            b_o[Tmod(maxindex+T-1)] = max(0, l[Tmod(maxindex+T-1)]-a_o[Tmod(maxindex+T-2)], r[maxindex]-a_o[maxindex])
            minmiles_odd = min(minmiles_odd, sum(a_o) + sum(b_o))
        end
        #println("Min miles = $minmiles_odd")
    elseif maxside == "L"
        for shift in 0:maxnum
            a_o[maxindex] = l[maxindex] - shift
            for t in [Tmod(t2) for t2 in maxindex+2:2:maxindex+T-1]
                b_o[Tmod(t-1)] = max(0, l[Tmod(t-1)]-a_o[Tmod(t-2)])
                a_o[t] = max(0, r[t]-b_o[Tmod(t-1)])
            end
            b_o[Tmod(maxindex+T-1)] = max(0, l[Tmod(maxindex+T-1)]-a_o[Tmod(maxindex+T-2)], r[maxindex]-a_o[maxindex])
            minmiles_odd = min(minmiles_odd, sum(a_o) + sum(b_o))
        end
        #println("Min miles = $minmiles_odd")
    end

    a_ev, b_ev = zeros(T), zeros(T)
    maxnum, maxindex, maxside = findmaxedgeweight_oddLs(r, l)
    minmiles_even = 1e10
    if maxside == "R"
        for shift in 0:maxnum
            b_ev[maxindex] = r[maxindex] - shift
            for t in [Tmod(t2) for t2 in maxindex+2:2:maxindex+T-1]
                a_ev[Tmod(t-1)] = max(0, r[Tmod(t-1)]-b_ev[Tmod(t-2)])
                b_ev[t] = max(0, l[t]-a_ev[Tmod(t-1)])
            end
            a_ev[Tmod(maxindex+T-1)] = max(0, r[Tmod(maxindex+T-1)]-b_ev[Tmod(maxindex+T-2)], l[maxindex]-b_ev[maxindex])
            minmiles_even = min(minmiles_even, sum(a_ev) + sum(b_ev))
        end
        #println("Min miles = $minmiles_even")
    elseif maxside == "L"
        for shift in 0:maxnum
            b_ev[maxindex] = l[maxindex] - shift
            for t in [Tmod(t2) for t2 in maxindex+2:2:maxindex+T-1]
                a_ev[Tmod(t-1)] = max(0, r[Tmod(t-1)]-b_ev[Tmod(t-2)])
                b_ev[t] = max(0, l[t]-a_ev[Tmod(t-1)])
            end
            a_ev[Tmod(maxindex+T-1)] = max(0, r[Tmod(maxindex+T-1)]-b_ev[Tmod(maxindex+T-2)], l[maxindex]-b_ev[maxindex])
            minmiles_even = min(minmiles_even, sum(a_ev) + sum(b_ev))
        end
        #println("Min miles = $minmiles_even")
    end

    a = a_o+a_ev
    b = b_o+b_ev

    return sum(a)+sum(b)
end
