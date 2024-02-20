
#Create output folders
if !(isdir("outputs"))
	mkdir("outputs")
	mkdir("outputs/static")
	#mkdir("outputs/static/figures")
elseif !(isdir("outputs/static/"))
	mkdir("outputs/static/")
	#mkdir("outputs/static/figures")
elseif !(isdir("outputs/static/figures"))
	mkdir("outputs/static/figures")
	#mkdir("outputs/static/figures")
end

if !(isdir(csvfoldername))
	mkdir(csvfoldername)
end
