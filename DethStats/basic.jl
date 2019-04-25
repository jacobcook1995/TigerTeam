#!/usr/bin/env julia
# basic.jl
# A script to do the basic parts of the data analysis, will most likely
# just be used for organisation
#
# Author: Jacob Cook
# Date: February 2019

using Distributions
using Plots
import PyPlot
using LsqFit

# Function to determine pearson smaple coefficent and associated p value
function pearson(xvec::Array{Float64,1},yvec::Array{Float64,1})
    # check two lists are the same length
    if length(xvec) != length(yvec)
        println("Error: Cannot find correlation of data sets that don't match.")
        error()
    end
    # need to include a step to remove NaNs from consideration
    inds = [] # find indices to delete
    for i = 1:length(xvec)
        if isnan(xvec[i]) || isnan(yvec[i])
            inds = vcat(inds,i)
        end
    end
    xvec = deleteat!(xvec,inds)
    yvec = deleteat!(yvec,inds)
    # Now calculate Pearson correlation coefficient
    xbarT = sum(xvec)/length(xvec)
    ybarT = sum(yvec)/length(yvec)
    a = 0
    b = 0
    c = 0
    for i = 1:length(xvec)
        a += (xvec[i] - xbarT)*(yvec[i] - ybarT)
        b += (xvec[i] - xbarT)^2
        c += (yvec[i] - ybarT)^2
    end
    r = a/sqrt(b*c)
    # Now calculate p value
    t = r*sqrt((length(xvec)-2)/(1-r^2)) # convert to t statistic
    d = TDist(length(xvec)-1) # make appropriate t distribution
    # Not 100% sure why the below is acceptable
    if t < 0.0
        P = 2*cdf(d,t) # use to calculate p value
    else
        P = 2*cdf(d,-t)
    end
    return(r,P)
end

# Make a function to do a linear least squares fit
function least2(xdata::Array{Float64,1},ydata::Array{Float64,1})
    # check two lists are the same length
    if length(xdata) != length(ydata)
        println("Error: Cannot find linear relation for data of different lengths")
        error()
    end
    # need to include a step to remove NaNs from consideration
    inds = [] # find indices to delete
    for i = 1:length(xdata)
        if isnan(xdata[i]) || isnan(ydata[i])
            inds = vcat(inds,i)
        end
    end
    xdata = deleteat!(xdata,inds)
    ydata = deleteat!(ydata,inds)
    # now define linear model
    @. model(x,p) = p[1] + p[2]*x
    # Need to remove NaN data
    p0 = [0.0,1.0]
    fit = curve_fit(model,xdata,ydata,p0)
    yint = coef(fit)[1]
    slope = coef(fit)[2]
    return(slope,yint)
end

# make a function to calculate standrad deviations
function standdev(data::Array{Float64,1},mean::Float64)
    n = length(data)
    a = 0
    for i = 1:n
        a += (data[i]-mean)^2
    end
    sd = sqrt(a/(n-1))
    return(sd)
end

function main()
    # Check there is a file of productions to be read
    infile = "../Input/soildata.csv"
    if ~isfile(infile)
        println("Error: No file of soil data to be read.")
        return(nothing)
    end
    # now read in 'Entropy productions'
    l = countlines(infile)
    w = 11
    soildata = Array{String,2}(undef,l,w)
    open(infile, "r") do in_file
        # Use a for loop to process the rows in the input file one-by-one
        k = 1
        for line in eachline(in_file)
            # parse line by finding commas
            L = length(line)
            comma = fill(0,w+1)
            j = 1
            for i = 1:L
                if line[i] == ','
                    j += 1
                    comma[j] = i
                end
            end
            comma[end] = L+1
            for i = 1:w
                soildata[k,i] = line[(comma[i]+1):(comma[i+1]-1)]
            end
            k += 1
        end
    end
    # Now a step to remove percentage signs using chop
    for i = 1:l
        # Chop everything but the NaNs
        if soildata[i,5] != "NaN"
            soildata[i,5] = chop(soildata[i,5])
        end
        if soildata[i,9] != "NaN"
            soildata[i,9] = chop(soildata[i,9])
        end
    end
    # Now have all the data in a string object
    idents = soildata[2:end,1] # seperate out identifying tags
    alt = parse.(Float64,soildata[2:end,2]) # Altitude in meters
    bd = parse.(Float64,soildata[2:end,3]) # Bulk density in g/cm^3
    erg = parse.(Float64,soildata[2:end,4]) # Ergostel (μg/g DW)
    loi = parse.(Float64,soildata[2:end,5]) # Loss on ignition percentage
    pH = parse.(Float64,soildata[2:end,6]) # pH
    TotN = parse.(Float64,soildata[2:end,7]) # percentage nitrogen
    TotC = parse.(Float64,soildata[2:end,8]) # percentage carbon total
    OrgC = parse.(Float64,soildata[2:end,9]) # percentage carbon Organic matter
    OrgNratio = parse.(Float64,soildata[2:end,10]) # Organic material to nitrogen ratio
    CNratio = parse.(Float64,soildata[2:end,11]) # Carbon nitrogen ratio
    # should find data ranges for grasslands/woodlands/grassslands-brigant
    fgrassrange = 1:49 # all grasslands
    woodrange = 50:94 # coal spoil (all) woodlands
    grassrange = 1:40 # coal spoil grasslands
    # Also make titles for output tables
    top = ["Alt m","BD","Ergosterol","LOI","pH","N%","C%","%C-OM","OM:N","C:N"]
    side = ["Altitude (m)","BD (g/cm^3)","Ergostel (ug/g DW)","LOI (%)","pH","Total N (%)","Total C (%)","%C (organic matter)","OM:N ratio","C:N ratio"]
    # First consider the woodland data
    wooddata = Array{String,2}(undef,11,11)
    wooddata[1,1] = "Woodland Soils"
    for i = 2:11
        wooddata[i,i] = "NaN"
    end
    wooddata[1,2:11] = top
    wooddata[2:11,1] = side
    # make new data set to do analysis on
    wood = hcat(alt[woodrange],bd[woodrange],erg[woodrange],loi[woodrange],pH[woodrange],TotN[woodrange],TotC[woodrange])
    wood = hcat(wood,OrgC[woodrange],OrgNratio[woodrange],CNratio[woodrange])
    # nested evaluation
    for i = 2:11
        for j = i+1:11
            r, P = pearson(wood[:,i-1],wood[:,j-1])
            wooddata[j,i] = "$r"
            wooddata[i,j] = "$P"
        end
    end
    outfile = "../Output/PearsonWood.csv"
    out_file = open(outfile, "w")
    for i = 1:11
        line = ""
        for j = 1:11
            line *= "$(wooddata[i,j]),"
        end
        line = line[1:end-1]
        line *= "\n"
        write(out_file,line)
    end
    close(out_file)
    # Then consider the full grassland data
    fgrassdata = Array{String,2}(undef,11,11)
    fgrassdata[1,1] = "All Grassland Soils"
    for i = 2:11
        fgrassdata[i,i] = "NaN"
    end
    fgrassdata[1,2:11] = top
    fgrassdata[2:11,1] = side
    # make new data set to do analysis on
    fgrass = hcat(alt[fgrassrange],bd[fgrassrange],erg[fgrassrange],loi[fgrassrange],pH[fgrassrange],TotN[fgrassrange],TotC[fgrassrange])
    fgrass = hcat(fgrass,OrgC[fgrassrange],OrgNratio[fgrassrange],CNratio[fgrassrange])
    # nested evaluation
    for i = 2:11
        for j = i+1:11
            r, P = pearson(fgrass[:,i-1],fgrass[:,j-1])
            fgrassdata[j,i] = "$r"
            fgrassdata[i,j] = "$P"
        end
    end
    outfile = "../Output/PearsonAllGrass.csv"
    out_file = open(outfile, "w")
    for i = 1:11
        line = ""
        for j = 1:11
            line *= "$(fgrassdata[i,j]),"
        end
        line = line[1:end-1]
        line *= "\n"
        write(out_file,line)
    end
    close(out_file)
    # Then consider just the spoil grassland data
    grassdata = Array{String,2}(undef,11,11)
    grassdata[1,1] = "Spoil Grassland Soils"
    for i = 2:11
        grassdata[i,i] = "NaN"
    end
    grassdata[1,2:11] = top
    grassdata[2:11,1] = side
    # make new data set to do analysis on
    grass = hcat(alt[grassrange],bd[grassrange],erg[grassrange],loi[grassrange],pH[grassrange],TotN[grassrange],TotC[grassrange])
    grass = hcat(grass,OrgC[grassrange],OrgNratio[grassrange],CNratio[grassrange])
    # nested evaluation
    for i = 2:11
        for j = i+1:11
            r, P = pearson(grass[:,i-1],grass[:,j-1])
            grassdata[j,i] = "$r"
            grassdata[i,j] = "$P"
        end
    end
    outfile = "../Output/PearsonSpoilGrass.csv"
    out_file = open(outfile, "w")
    for i = 1:11
        line = ""
        for j = 1:11
            line *= "$(grassdata[i,j]),"
        end
        line = line[1:end-1]
        line *= "\n"
        write(out_file,line)
    end
    close(out_file)
    # Now that the data has been output there are two graphs that seem worth plotting
    pyplot(dpi=150)
    xplotrange = 0.65:0.05:1.7
    scatter([bd[grassrange]],[erg[grassrange]],legend=false,xlims=(0.6,1.8),ylims=(0.6,14.0))
    plot!(xlabel="Bulk Density (g.cm^-3)",ylabel="Ergosterol (ug.g^-1)",xticks=0.6:0.2:1.8,yticks=0.0:2.0:14.0)
    slope, yint =  least2(bd[grassrange],erg[grassrange])
    plot!(xplotrange,yint.+slope.*xplotrange,color=1)
    savefig("../Output/ErgosterolvsBulkDensity.png")
    xplotrange = 125.0:25.0:375.0
    scatter([alt[grassrange]],[bd[grassrange]],legend=false,xlims=(100.0,400.0),ylims=(0.6,1.8))
    plot!(xlabel="Quadrat Altitude (m)",ylabel="Bulk Density (g.cm^-3)",xticks=100.0:50.0:400.0,yticks=0.6:0.2:1.8)
    slope, yint =  least2(alt[grassrange],bd[grassrange])
    plot!(xplotrange,yint.+slope.*xplotrange,color=1)
    savefig("../Output/BulkDensityvsAltitude.png")
    # Now to make the bar charts for grassland data
    BD = bd[fgrassrange]
    ERG = erg[fgrassrange]
    labels = idents[fgrassrange]
    tokens = Array{String,1}(undef,length(labels))
    # Reduce these labels to two letter idetifiers
    for i = 1:length(labels)
        token = labels[i]
        tokens[i] = token[1:2]
    end
    # New data structures for each site
    sites = unique(tokens) # slow but okay method of obtaining unqiue entries
    bulk = zeros(length(sites))
    ergo = zeros(length(sites))
    erbulk = zeros(length(sites))
    erergo = zeros(length(sites))
    # now use for loop to find ranges for each site
    ranges = zeros(Int64,length(sites),2)
    count = 0
    for i = 1:length(sites)
        all = false
        ranges[i,1] = count+1
        token = tokens[count+1]
        while all == false
            count += 1
            if count+1 > length(labels) || tokens[count+1] != token
                all = true
            end
        end
        ranges[i,2] = count
    end
    # Now use a loop to calculate means and standard deviations
    for i = 1:length(sites)
        # filter relevant data
        filterg = filter(y->!isnan(y),ERG[ranges[i,1]:ranges[i,2]])
        filtbd = filter(y->!isnan(y),BD[ranges[i,1]:ranges[i,2]])
        # find means
        ergo[i] = sum(filterg)/length(filterg)
        bulk[i] = sum(filtbd)/length(filtbd)
        # use my function to find standard deviations
        erergo[i] = standdev(filterg,ergo[i])
        erbulk[i] = standdev(filtbd,bulk[i])
    end
    # Now plot two bar charts
    bar(sites,ergo,yerr=erergo,legend=false,ylabel="Ergosterol (ug.g^-1)",ylims=(0,18),yticks=0.0:2.0:18.0,markercolor=:black,markerstrokecolor=:black)
    savefig("../Output/ErgosterolBar.png")
    bar(sites,bulk,yerr=erbulk,legend=false,ylabel="Bulk Denisty (g.cm^-1)",ylims=(0,1.8),yticks=0.0:0.2:1.8,markercolor=:black,markerstrokecolor=:black)
    savefig("../Output/BulkDenistyBar.png")
    # Now should count number of unique tokens
    return(nothing)
end

@time main()
