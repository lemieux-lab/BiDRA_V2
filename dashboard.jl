using GenieFramework
using CSV, DataFrames
include("juBiDRA.jl")
include("juPlot.jl")
@genietools

Genie.config.cors_headers["Access-Control-Allow-Origin"]  =  "*"
Genie.config.cors_headers["Access-Control-Allow-Headers"] = "Content-Type"
Genie.config.cors_headers["Access-Control-Allow-Methods"] = "GET,POST,PUT,DELETE,OPTIONS"
Genie.config.cors_allowed_origins = ["*"]

const FILE_PATH = "upload/" ## Where to save user dataset (CSV)
mkpath(FILE_PATH)

const ANALYSIS_PATH = "analysis/"
mkpath(ANALYSIS_PATH)

@in responseSelected = "Ascending" ## Default selection of response tye
@out title = "BiDRA Dashboard" ## Main title on the page
@out responseOption = ["Ascending", "Descending"] ## Options of response typ for the user

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
  responseType = Genie.Requests.postpayload(:resp)
  
  #uniqueID = "0203451612" ## For development purposes
  #responseType = "Descending" ## For development purposes
  println("##### $uniqueID Analysis #####")
  global dataset_path = "$ANALYSIS_PATH$uniqueID/"
  mkpath(dataset_path)

  chrono = @elapsed dataset = getDataset(uniqueID)
  datasetExp = groupby(dataset, :id)
  println("--- Dataset uploaded ($chrono sec.)")

  expId = unique(dataset.id)
  N = length(expId)
  println("--- There are $N individual analysis to analyse")

  println("--- Starting inference ")
  for i in 1:N
    subsetId = expId[i];
    println("------ ($i) $subsetId")

    subsetData = DataFrame(datasetExp[i]);
    chrono = @elapsed posterior = doInference(subsetData, subsetId, responseType);
    println("--------- Posterior Inferred ($chrono sec.)")


    chrono = @elapsed plotDoseResponse(subsetData, posterior, subsetId);
    println("--------- Inference Plotted ($chrono sec.)")

    
  end

  println("--- DONE")
  
  #"Hello $(Genie.Requests.postpayload(:resp, "Anon"))"
end



### Views renderer
@page("/", "views/app.jl.html")
Server.isrunning() || Server.up()
