
using Plots, CSV, DataFrames

filename = "outputs/driverequity/driverhours_combined.csv"

data = CSV.read(filename, DataFrame)

noequitydata = filter(row -> row.formulation == "homogeneousdeadlines", data)
equitydata = filter(row -> row.formulation == "heterogeneous", data)

p1 = histogram(noequitydata[:, :hours], bins = 0:2:40, ylims=(0, 600),  xlabel="X", ylabel="Driver count") 
savefig(p1, "noequity.png")

p2 = histogram(equitydata[:, :hours], bins = 0:2:40, ylims=(0, 600)) 
savefig(p2, "equity.png")