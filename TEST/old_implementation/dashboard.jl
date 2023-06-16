using GenieFramework, HTTP
using CSV, DataFrames
@genietools

include("juBiDRA.jl") ## utils functions
include("juPlot.jl") ## plotting functions

Genie.config.cors_headers["Access-Control-Allow-Origin"]  =  "*"
Genie.config.cors_headers["Access-Control-Allow-Headers"] = "Content-Type"
Genie.config.cors_headers["Access-Control-Allow-Methods"] = "GET,POST,PUT,DELETE,OPTIONS"
Genie.config.cors_allowed_origins = ["*"]

const FILE_PATH = "upload/" ## Where to save user dataset (CSV)
mkpath(FILE_PATH)

const ANALYSIS_PATH = "analysis/" ## Where analysis results are saved
mkpath(ANALYSIS_PATH)

@in responseSelected = "Ascending" ## Default selection of response tye
@out responseOption = ["Ascending", "Descending"] ## Options of response typ for the user
@out isprocess = false
global uniqueID = "" ## Initiate


## When a new dataset is submitted
route("/upload", method = POST) do
  files = Genie.Requests.filespayload()

  for f in files ## eventually could add mutiple files options
    global uniqueID = getUniqueID()
    global experimentName = split(f[2].name, ".")[1]
    write(joinpath(FILE_PATH, "$uniqueID.csv"), f[2].data)
  end

  if length(files) == 0
      @info "No file uploaded"
  end
  return "Upload finished"
end

## When a new analysis is requested
route("/", method = POST) do
  println(getUniqueID)
  responseType_in = Genie.Requests.postpayload(:resp)
  ## 0 : Descending 
  ## 1 : Ascending 
  responseType = responseType_in == "Ascending" ? 1 : 0
  
  #uniqueID = "0203451612" ## For development purposes
  #responseType = "Descending" ## For development purposes
  println("##### $uniqueID Analysis #####")
  global dataset_path = "$ANALYSIS_PATH$uniqueID/"
  mkpath(dataset_path)
  println("$ANALYSIS_PATH$uniqueID/")

  chrono = @elapsed dataset = getDataset(uniqueID)
  #datasetExp = groupby(dataset, :id)
  println("--- Dataset uploaded ($chrono sec.)")

  expId = unique(dataset.id)
  N = length(expId)
  println("--- There are $N individual analysis to analyse")

  SD_ALL = []
  ΔWAICₖ_ALL = []
  posterior_ALL = DataFrame()

  println("--- Starting inference ")
  for i in 1:N
    subsetId = expId[i];
    println("------ ($i) $subsetId")

    subsetDataset = filter(:id => x -> x == subsetId, dataset)
    chrono = @elapsed posterior = doInference(subsetData, subsetId, responseType);
    posterior[!, :exp_id] = repeat([subsetId], size(posterior)[1])
    posterior_ALL =  vcat(posterior_ALL, posterior)
    println("--------- Posterior Inferred ($chrono sec.)")

    chrono = @elapsed plotDoseResponse(subsetData, posterior, subsetId, responseType);
    println("--------- Inference Plotted ($chrono sec.)")

    waicₖ_bidra, waicₖ_line = get_waicₖ(subsetData, posterior, responseType)
    ΔWAICₖ = waicₖ_bidra - waicₖ_line
    ΔWAICₖ_ALL = vcat(ΔWAICₖ_ALL, ΔWAICₖ)
    SD_ALL = vcat(SD_ALL, std(subsetData.response))
    println("--------- Metrics calculated")
 
  end

  dataset_metrics = DataFrame(exp_id=expId, SD=SD_ALL, ΔWAICₖ=ΔWAICₖ_ALL)
  dataset_metrics[!, :completeness] = map(r -> r ≥ 20 ? 1 : 2, dataset_metrics.SD)
  dataset_metrics[!, :WAICₖ_model] = map(r -> r < 0 ? 1 : 2, dataset_metrics.ΔWAICₖ)
  dataset_metrics[!, :group] = dataset_metrics.WAICₖ_model .+ dataset_metrics.completeness .- 1

  sort!(dataset_metrics, :group, rev=true)

  chrono = @elapsed plotGroupAssignation(dataset_metrics);
  println("--------- Group assignation (or posterior informative potentiel) plotted ($chrono sec.)")

  if N ≥ 3
    ## ranking
  elseif  N == 2
    ## comparison
    chrono = @elapsed plotPairedComparison(posterior_ALL, expId)
    println("--------- Pairs of posterior comparison plotted ($chrono sec.)")
  end

  println("--- DONE")

  
  ## Return file download
  getPwd = pwd()
  cd(ANALYSIS_PATH)
  
  chrono = @elapsed fn = createZip(uniqueID, getPwd)
  println("--- Analysis ZIP created ($chrono sec.)")

  @out isprocess = false

  HTTP.Response(200, ["Content-Type" => "application/zip"], body = read(ANALYSIS_PATH*fn))
end


### Views renderer
@page("/", "views/app.jl.html")

## Return true or launch server 
Server.isrunning() || Server.up()
