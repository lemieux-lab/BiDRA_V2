using Turing, MCMCChains
using Statistics, StatsBase
using DataFrames, CSV
using CairoMakie

include("juPlot.jl")

function getUniqueID()
    return Dates.format(Dates.now(), "SSmmMMddHH")
end

function replaceInDf(df::DataFrame, col::Symbol, from::String, to::String)
    df[!, col] .= replace.(df[:, col], to => from)
    return df
end

function readCSV(fn::String)
    data = CSV.read(fn, header=false, DataFrame)
    return data
end

function getDataset(uniqueID::String)
    path = "$FILE_PATH/$uniqueID.csv"
    dataset = readCSV(path)
    rename!(dataset, [:concentration, :response, :id])

    dataset[!, :id] = string.(dataset.id)

    data_df = replaceInDf(dataset, :id, " ",  "_")
    data_df = replaceInDf(dataset, :id, "/",  "_")
    data_df = replaceInDf(dataset, :id, ",",  "_")
    data_df = replaceInDf(dataset, :id, "(",  "_")
    data_df = replaceInDf(dataset, :id, ")",  "_")

    return data_df
end

function llogistic(param::Array) 
    HDR, LDR, ic50, slope = param 

    return x -> HDR + ((LDR - HDR) / (1 + 10^(slope * (x - ic50))))
end

@model function sigmoid_BIDRA(xs::Array, ys::Array, HDRmixture::MixtureModel, LDRμ::Int) 
    HDR ~ HDRmixture
    LDR ~ Normal(LDRμ, 10)
    ic50 ~ Normal(0,10)
    slope ~ LogNormal(0.5, 1)
    
    σ ~ LogNormal(1, 1)

    for i in 1:length(xs)
        f = llogistic([HDR, LDR, ic50, slope])
        ys[i] ~ Normal(f(xs[i]), σ)
    end
end

@model function continuous_BIDRA(xs::Array, ys::Array, LDRμ::Int)
    line ~ Normal(LDRμ, 10)
    σ ~ LogNormal(1, 1)

    for i in 1:length(xs)
        ys[i] ~ Normal(line, σ)
    end
end

function doInference(expData::DataFrame, id::String, responseType::String)
    local HDRmixture
    local LDRμ

    λ = [0.4, 0.5, 0.1]
    if responseType == "Ascending"
        HDRmixture = MixtureModel([SkewNormal(100, 10, 1), Uniform(0, 100), SkewNormal(0, 10, -1)], λ)
        LDRμ = 0
    elseif responseType == "Descending"
        HDRmixture= MixtureModel([SkewNormal(0, 10, 1), Uniform(0, 100), SkewNormal(100, 20, -5)], λ)
        LDRμ = 100
    end

    nChain = 4
    nIte = 1000
    nAdapt = 1000
    δ = 0.65
    Turing.setprogress!(false)

    modelₛ = sigmoid_BIDRA(expData.concentration, expData.response, HDRmixture, LDRμ)
    modelₗ = continuous_BIDRA(expData.concentration, expData.response, LDRμ)
    sampler = NUTS(nAdapt, δ)

    timeₛ = @elapsed chainsₛ = sample(modelₛ, sampler, MCMCThreads(), nIte, nChain)
    timeₗ = @elapsed chainsₗ = sample(modelₗ, sampler, MCMCThreads(), nIte, nChain)

    chainsCombined = DataFrame(chainsₛ)[:, [:HDR, :LDR, :ic50, :slope]]
    chainsCombined[!, :σₛ] = vcat(chainsₛ[:σ]...)
    chainsCombined[!, :cont] = vcat(chainsₗ[:line]...)
    chainsCombined[!, :σₗ] = vcat(chainsₗ[:σ]...)

    ## save results
    posterior_path = dataset_path*"posterior/"
    mkpath(posterior_path)

    fileName = posterior_path*"posterior_expID_$id.csv"
    CSV.write(fileName, chainsCombined)

    return chainsCombined
end


function download(filepath::String; root)::HTTP.Response 
    return serve_static_file(filepath; root=root, download=true) 
end 

