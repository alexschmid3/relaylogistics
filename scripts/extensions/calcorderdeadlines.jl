
function calcorderdeadlines(shortesttriptimes)

    orderdeadline = Dict()

    for i in orders
        timetodeadline = tstep * ceil(shortesttriptimes[i]/tstep) * deadlineasmultipleofshortestpath
        orderdeadline[i] = orderOriginalStartTime[i] + timetodeadline
    end

    return orderdeadline

end
