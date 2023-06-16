using CSV, DataFrames
using Turing, MCMCChains

include("juBiDRA.jl")

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

function resetOuput()
  return "", "", "", "", "", ""
end

function getResponseTypeInt(inVal)
    return inVal == "Ascending" ? 1 : 0
end

function replaceInDf(df::DataFrame, col::Symbol, from::String, to::String)
  df[!, col] .= replace.(df[:, col], to => from)
  return df
end

function readCSV(fn::String)
  data = CSV.read(fn, header=false, DataFrame)
  return data
end

function getDataset(filePath::String, uniqueID::String)
  path = "$filePath/$uniqueID.csv"
  dataset = readCSV(path)
  rename!(dataset, [:concentration, :response, :id])

  dataset[!, :id] = string.(dataset.id)

  data_df = replaceInDf(dataset, :id, " ",  "_")
  data_df = replaceInDf(dataset, :id, "/",  "_")
  data_df = replaceInDf(dataset, :id, ",",  "_")
  data_df = replaceInDf(dataset, :id, "(",  "_")
  data_df = replaceInDf(dataset, :id, ")",  "_")

  expId = unique(dataset.id)
  N = length(expId)

  return data_df, expId, N
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

  model_toRun_bidra = model_BIDRA(expData.concentration, expData.response, HDRmixture, LDRμ, responseType)
  model_toRun_line = model_line(expData.concentration, expData.response, LDRμ)
  sampler = NUTS(nAdapt, δ)

  chains_bidra = DataFrame(sample(model_toRun_bidra, sampler, MCMCThreads(), nIte, nChain))[:, [:HDR, :LDR, :ic50, :slope, :σ]]
  chains_line = DataFrame(sample(model_toRun_line, sampler, MCMCThreads(), nIte, nChain))[:, [:line, :σ]]

  chainsCombined = copy(chains_bidra)
  chainsCombined[!, :line] = vcat(chains_line[:, :line]...)
  chainsCombined[!, :σ_line] = vcat(chains_line[:, :σ]...)

  ## save results 
  fileName = posterior_path*"posterior_expID_$id.csv"
  CSV.write(fileName, chains_bidra)

  return chainsCombined
end

function get_bidra_μ(posterior_df::DataFrame, experimental_dose::Vector{Float64}, responseType::Int)
  param_posterior = posterior_df[:, [:HDR, :LDR, :ic50, :slope, :σ]]
  f_ite = map(row -> llogistic(Array(row), responseType), eachrow(param_posterior))
  μ =  mapreduce(f -> f.(experimental_dose), hcat, f_ite) ## μ of likelihood function for each dose at each Bayeisna iteration

  return μ
end

function get_log_pdf_bidra(data_df::DataFrame, posterior::DataFrame, responseType::Int, logTrans::Bool)
  experimental_dose = data_df.concentration
  experimental_response = data_df.response
  n = nrow(posterior)

  μ = get_bidra_μ(posterior, experimental_dose, responseType)
  σ = posterior[:, :σ]

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
  σ = posterior[:, :σ_line]

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

function getSortedRankedMatrix(values_matrix::Matrix{Float64}, N::Int, n_samples::Int, sorted_expId::Vector, paramRev::Bool)
  ranks_matrix = mapslices(x -> ordinalrank(x, rev=!paramRev), values_matrix, dims=2)
  counts_rank = countmap.(eachcol(ranks_matrix))

  ## Calculate prob for each rank (rank x compound)
  counts_matrix = zeros(Float64, N, N)
      
  for idx in 1:N
    exp_count = counts_rank[idx]

    for (rank, count) in pairs(exp_count)
      counts_matrix[rank, idx] = count/n_samples
    end
  end
  return counts_matrix[:, sorted_expId]
end

function rankingAnalysis(posterior::DataFrame, N::Int, expId::Vector, paramRank::Int, responsetype::Int)
  if paramRank == 0
    param = :ic50
    paramRev = true ## Rank1 = smallest
  elseif paramRank == 1
    param = :HDR
    if responseType == 0
      paramRev = true ## Rank1 = smallest
    elseif  responseType == 1
      paramRev = false ## Rank1 = highest
    end
  end

  param_posterior = posterior[:, [param, :exp_id]]
  n_samples = Int(size(param_posterior)[1] / N)

  ## rank by median value (to get some sort of diagonal)
  param_median = combine(groupby(param_posterior, :exp_id), param => median => :median)
  sort!(param_median, :median, rev=paramRev) 
  sorted_pos = [findfirst(x -> x == e, expId) for e in param_median.exp_id]

  ## Create matrices (sampling x compound)
  values_matrix = reshape(param_posterior[:, param], (n_samples, N))
  sorted_rank_matrix = getSortedRankedMatrix(values_matrix, N, n_samples, sorted_pos, paramRev)

  plotRanking(sorted_rank_matrix, posterior, sorted_pos, expId, N, n_samples, paramRank)
end

function createZip(uniqueID::String, resetPwd::String)
  writerName = "results_$uniqueID.zip"
  cmmd = `zip -r $writerName $uniqueID/`
  run(cmmd)

  cd(resetPwd)
  rm("$UPLOAD_PATH/$uniqueID.csv")
  println("--- Uploaded dataset deleted from server")

  rm("$ANALYSIS_PATH/", recursive=true)
  println("--- Analysis results dir deleted from server")
  return writerName
end