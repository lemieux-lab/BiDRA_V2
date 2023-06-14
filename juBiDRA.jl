using Turing, MCMCChains
using Statistics, StatsBase
using DataFrames, CSV
using CairoMakie
using ZipFile

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

function llogistic(param::Array, responseType::Int) 
    HDR, LDR, ic50, slope = param 

    if responseType == 1
        return x -> LDR + ((HDR - LDR) / (1 + 10^(slope * (ic50 - x))))
    else
        return x -> HDR + ((LDR - HDR) / (1 + 10^(slope * (x - ic50))))
    end
end

@model function model_BIDRA(xs::Array, ys::Array, HDRmixture::MixtureModel, LDRμ::Int, responseType::Int) 
    HDR ~ HDRmixture
    LDR ~ Normal(LDRμ, 10)
    ic50 ~ Normal(0,10)
    slope ~ LogNormal(0.5, 1)
    
    σ ~ LogNormal(1, 1)

    for i in 1:length(xs)
        f = llogistic([HDR, LDR, ic50, slope], responseType)
        ys[i] ~ Normal(f(xs[i]), σ)
    end
end

@model function model_line(xs::Array, ys::Array, LDRμ::Int)
    line ~ Normal(LDRμ, 10)
    σ ~ LogNormal(1, 1)

    for i in 1:length(xs)
        ys[i] ~ Normal(line, σ)
    end
end

function doInference(expData::DataFrame, id::String, responseType::Int)
    local HDRmixture
    local LDRμ

    λ = [0.4, 0.5, 0.1]
    if responseType == 1
        HDRmixture = MixtureModel([SkewNormal(100, 10, 1), Uniform(0, 100), SkewNormal(0, 10, -1)], λ)
        LDRμ = 0
    else
        HDRmixture= MixtureModel([SkewNormal(0, 10, 1), Uniform(0, 100), SkewNormal(100, 20, -5)], λ)
        LDRμ = 100
    end

    nChain = 4
    nIte = 1000
    nAdapt = 1000
    δ = 0.65
    Turing.setprogress!(false)

    modelᵦ = model_BIDRA(expData.concentration, expData.response, HDRmixture, LDRμ, responseType)
    modelₗ = model_line(expData.concentration, expData.response, LDRμ)
    sampler = NUTS(nAdapt, δ)

    timeₛ = @elapsed chainsᵦ = sample(modelᵦ, sampler, MCMCThreads(), nIte, nChain)
    timeₗ = @elapsed chainsₗ = sample(modelₗ, sampler, MCMCThreads(), nIte, nChain)

    chainsCombined = DataFrame(chainsᵦ)[:, [:HDR, :LDR, :ic50, :slope]]
    chainsCombined[!, :σᵦ] = vcat(chainsᵦ[:σ]...)
    chainsCombined[!, :line] = vcat(chainsₗ[:line]...)
    chainsCombined[!, :σₗ] = vcat(chainsₗ[:σ]...)

    ## save results
    posterior_path = dataset_path*"posterior/"
    mkpath(posterior_path)

    fileName = posterior_path*"posterior_expID_$id.csv"
    CSV.write(fileName, chainsCombined)

    return chainsCombined
end

function get_bidra_μ(posterior_df::DataFrame, experimental_dose::Vector{Float64}, responseType::Int)
    param_posterior = posterior_df[:, [:HDR, :LDR, :ic50, :slope, :σᵦ]]
    f_ite = map(row -> llogistic(Array(row), responseType), eachrow(param_posterior))
    μ =  mapreduce(f -> f.(experimental_dose), hcat, f_ite) ## μ of likelihood function for each dose at each Bayeisna iteration

    return μ
end

function get_log_pdf_bidra(data_df::DataFrame, posterior::DataFrame, responseType::Int, logTrans::Bool)
    experimental_dose = data_df.concentration
    experimental_response = data_df.response
    n = nrow(posterior)

    μ = get_bidra_μ(posterior, experimental_dose, responseType)
    σ = posterior[:, :σᵦ]

    if logTrans
        return mapreduce(i -> vcat(logpdf.(Normal.(μ[:, i], σ[i]), experimental_response)...), hcat, 1:n)
    else
        return mapreduce(i -> vcat(pdf.(Normal.(μ[:, i], σ[i]), experimental_response)...), hcat, 1:n)
    end
end

function get_log_pdf_line(data_df::DataFrame, posterior::DataFrame, logTrans::Bool)
    experimental_response = data_df.response
    n = nrow(posterior)

    μ = posterior[:, :line]
    σ = posterior[:, :σₗ]

    if logTrans
        return mapreduce(i -> vcat(logpdf.(Normal(μ[i], σ[i]), experimental_response)...), hcat, 1:n)
    else
        return mapreduce(i -> vcat(pdf.(Normal(μ[i], σ[i]), experimental_response)...), hcat, 1:n)
    end
end

function get_waicₖ(data_df::DataFrame, posterior::DataFrame, responseType::Int)
    lpd_bidra = get_log_pdf_bidra(data_df, posterior, responseType, false)
    lpd_line = get_log_pdf_line(data_df, posterior, false)

    lppd_bidra = sum(log.(sum.(eachrow(lpd_bidra)) ./ size(lpd_bidra)[2]))
    lppd_line = sum(log.(sum.(eachrow(lpd_line)) ./ size(lpd_line)[2]))

    waicₖ_bidra = -2 * (lppd_bidra - 5)
    waicₖ_line = -2 * (lppd_line - 2)

    return waicₖ_bidra, waicₖ_line
end

function createZip(uniqueID::String, resetPwd::String)
    writerName = "results_$uniqueID.zip"
    cmmd = `zip -r $writerName $uniqueID/`
    run(cmmd)

    cd(resetPwd)
    rm(FILE_PATH*"$uniqueID.csv")
    println("--- Uploaded dataset deleted from server")

    rm(dataset_path, recursive=true)
    println("--- Analysis results dir deleted from server")
    return writerName
end 

