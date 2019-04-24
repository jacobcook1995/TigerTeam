#!/usr/bin/env julia
# basic.jl
# A script to do the basic parts of the data analysis, will most likely
# just be used for organisation
#
# Author: Jacob Cook
# Date: February 2019

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
    OrgNratio = parse.(Float64,soildata[2:end,10])
    CNratio = parse.(Float64,soildata[2:end,11])
    # Data is now appropriately processed
    
    return(nothing)
end

@time main()
