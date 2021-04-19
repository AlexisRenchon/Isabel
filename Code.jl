# First, make sure your working directory is set

# Load packages
# CSV is a package to read and write text files
# DataFrames is a package to work with tabular data 
# Dates is a package to work with Dates 
using CSV, DataFrames, Dates

# This will just read your input data
data = CSV.read(joinpath("Input","projectDAMM.csv"), DataFrame)

# This creates a vector of datetime from your data. 
# Note the . (DateTime.()), this is to use the function on all your elements
datetime = DateTime.(data.YR, data.MO, data.DAY, data.HR, data.MIN, data.SEC)

# This load specific columns of your data
SWC = data."SWC (%)"

# Here I replace your NaN by 0, I get rid of the 0 after
SWC = replace(SWC, NaN=>0)
SWC = SWC./100 # convert to m³ m⁻³
SWC_f = SWC[SWC .> 0] # delete the 0, this is to only keep time where you have data

Tsoil = data."T(K)" .- 273.15 # convert to °C, note again the . (.-) 
Tsoil_f = Tsoil[SWC .> 0] # only keep data when you have SWC data

Rsoil =  data."F_CO2 (u mol)"
Rsoil_f = Dep_var = Rsoil[SWC .> 0]

datetime_f = datetime[SWC .> 0]

# This creates a matrix of 2 columns
Ind_var = hcat(Tsoil_f, SWC_f)

# DAMM model package
using DAMMmodel

# Test with some default parameters
p = Param_ini = [1e8, 62, 3.46e-8, 2.0e-3, 0.34, 0.0125] 

# This estimates Rsoil with your Tsoil and SWC with default param
output1 = DAMM(Ind_var, p)

# You want to fit parameters to your Rsoil data
# LsqFit will find the parameters that minimize the error between data and model
using LsqFit

# This fit the DAMM parameters to your data, it minimize the square of the errors between modeled Rsoil and observed Rsoil
fit = curve_fit(DAMM, Ind_var, Dep_var, Param_ini)
Param_fit = coef(fit) # these are the new values of your parameters
Modeled_data = DAMM(Ind_var,Param_fit) # that's the modeled Rsoil estimated for your Tsoil and SWC data, you could also use it on other Tsoil and SWC

# This is a basic plotting package
using Plots

scatter(datetime_f, Modeled_data) 
scatter!(datetime_f, Rsoil_f) # the ! (scatter!()) is to add a plot on the previous one

# This creates a DataFrame with your data and the output
data_out = DataFrame(datetime = datetime_f,
		     SWC = SWC_f,
		     Tsoil = Tsoil_f,
		     Rsoil = Rsoil_f,
		     RsoilDAMM = Modeled_data)

# This save an output file 
CSV.write(joinpath("Output","output.csv"), data_out)

# Read met data file
met = CSV.read(joinpath("Input","Met_Julia_03122021.csv"), DataFrame;
	       # missingstrings = ["NAN", ""],
	       header = 2,
	       datarow = 7079, # when you start collecting SWC data
	       dateformat = "mm/dd/yyyy HH:MM")

# Calculate fluxes with DAMM fitted params and met
Ind_var_met = hcat(met."TC_Avg(1)", met.VWC_Avg) # need to get rid of missing
Param_fit[5] = 0.56
Model_Rsoil_met = DAMM(Ind_var_met, Param_fit)

using GLMakie
using UnicodeFun
using PlotUtils: optimize_ticks
Dtime_all_rata = datetime2rata.(met.TIMESTAMP)
dateticks = optimize_ticks(met.TIMESTAMP[1], met.TIMESTAMP[end])[1]

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
ax3 = Axis(fig[3,1])

p1 = GLMakie.scatter!(ax1, Dtime_all_rata, Model_Rsoil_met, color = :black) # there are periods with VWC = 0 which is impossible
p2 = GLMakie.scatter!(ax2, Dtime_all_rata, met."TC_Avg(1)", color = :red)
p3 = GLMakie.scatter!(ax3, Dtime_all_rata, met.VWC_Avg, color = :blue)

# Attributes of Axis, e.g., ax1.attributes 
ax1.ylabel = to_latex("R_{soil} (\\mumol m^{-2} s^{-1})");
ax1.xticks[] = (datetime2rata.(dateticks) , Dates.format.(dateticks, "yyyy"));
ax1.xticklabelsvisible = false
ax2.ylabel = to_latex("T_{soil} (°C)");
ax2.xticks[] = (datetime2rata.(dateticks) , Dates.format.(dateticks, "yyyy"));
ax2.xticklabelsvisible = false
ax3.ylabel = to_latex("\\theta (m^3 m^{-3})");
ax3.xticks[] = (datetime2rata.(dateticks) , Dates.format.(dateticks, "yyyy"));

# Attributes of plots, e.g., p1.attributes
p1.strokewidth = 0
p2.strokewidth = 0
p3.strokewidth = 0

fig

# 3D plot
fig2 = Figure()
ax3D = Axis3(fig2[1,1])
p3D = GLMakie.scatter!(ax3D, Tsoil_f, SWC_f, Rsoil_f, markersize = 5000, color = :black)

L = 40 # resolution
x = collect(range(1, length=L, stop=1))
[append!(x, collect(range(i, length=L, stop=i))) for i = 2:40]
x = reduce(vcat,x)
y = collect(range(0.01, length=L, stop=0.40))
y = repeat(y, outer=L)
x_range = hcat(x, y)

x = Int.(x_range[:, 1])
y_ax = collect(range(0.01, length=L, stop=0.40))
y = collect(range(1, length=L, stop=L))
y = repeat(y, outer=L)
y = Int.(y)
x_ax = collect(range(1, length=L, stop=L))

using SparseArrays
surface!(ax3D, x_ax, y_ax, Matrix(sparse(x, y, DAMM(x_range, Param_fit))), colormap = Reverse(:Spectral), transparency = true, alpha = 0.2, shading = false)

wireframe!(ax3D, x_ax, y_ax, Matrix(sparse(x, y, DAMM(x_range, Param_fit))), overdraw = true, transparency = true, color = (:black, 0.1));


ax3D.xlabel = to_latex("T_{soil} (°C)");
ax3D.ylabel = to_latex("\\theta (m^3 m^{-3})");
ax3D.zlabel = to_latex("R_{soil} (\\mumol m^{-2} s^{-1})");


