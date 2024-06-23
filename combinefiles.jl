
using CSV, DataFrames

figuredirectory = string("outputs/online/stability/")

#Get list of files
filelist = readdir(figuredirectory)

#Get the file types to combine
filefirstwords = [file[1:2] for file in filelist]
#filefirstwords = [file[1:findfirst("_", file)[1]-1] for file in filelist]
filetypes = unique!(filefirstwords)
dataDict = Dict()
for t in filetypes
	dataDict[t] = []
end
for file in filelist
	filetype = file[1:2]
	push!(dataDict[filetype], file)
end

#Combine each file type to one csv
for t in filetypes
	
	#Name of combo file
	combofile = string(figuredirectory, "/", t,"_combined.csv")

	#Combine files and save to combofile
	for file in dataDict[t]
		filenum = file[4:findfirst("_", file)[1]-1]
		#Add all rows (including header) from first file
		if file == dataDict[t][1]
			open(string(figuredirectory, "/", file)) do input
			    open(combofile, "a") do output
			        for line in eachline(input)
			            println(output, string(line))
			        end
			    end
			end
		#Drop header row from other files
		else
			open(string(figuredirectory, "/", file)) do input
			    open(combofile, "a") do output
			        for line in Iterators.drop(eachline(input), 1)
			            println(output, string(line))
			        end
			    end
			end
		end
	end

end 
