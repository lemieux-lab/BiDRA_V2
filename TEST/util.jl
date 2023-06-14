function create_storage_dir(name)
    try
      mkdir(joinpath(@__DIR__, name)) 
    catch 
      @warn "directory already exists" 
    end
    return joinpath(@__DIR__, name)
end

function getUniqueID()
    return Dates.format(Dates.now(), "SSmmMMddHH")
end

function getResponseTypeInt(inVal)
    return inVal == "Ascending" ? 1 : 0
end