#!/usr/bin/env julia
# basic.jl
# A script to do the basic parts of the data analysis, will most likely
# just be used for organisation
#
# Author: Jacob Cook
# Date: February 2019

using Distributions

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
    return(nothing)
end

@time main()
