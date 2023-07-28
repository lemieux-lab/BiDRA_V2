using GenieFramework
using HTTP
@genietools

include("util.jl")
include("juPlot.jl")

@appname BiDRA_V2

Genie.config.cors_headers["Access-Control-Allow-Origin"]  =  "*"
Genie.config.cors_headers["Access-Control-Allow-Headers"] = "Content-Type"
Genie.config.cors_headers["Access-Control-Allow-Methods"] = "GET,POST,PUT,DELETE,OPTIONS"
Genie.config.cors_allowed_origins = ["*"]

const UPLOAD_PATH = create_storage_dir("Backend_Upload")
const RESULTS_PATH = create_storage_dir("Analysis_Results")

#GenieSession.__init__()

@app begin
  @in process = false
  @out inputID = getUniqueID()

  @in select = "Ascending"
  @in option = ["Ascending", "Descending"]

  @out output1 = ""
  @out output2 = ""

  @onbutton process begin

    load = true

    output1 = select 
    output2 = "Your DR experiments are beeing analyzed. Once the analysis process is completed, the results will be automatically downloaded. Thank you for your patience!"
  end
end

function ui(uniqueID)
  row(cell(class = "st-module", [
    uploader(name="fileUpload", label="Upload Dataset", accept=".csv", multiple=false, method="POST", url="/upload/$uniqueID", autoupload=true)

    StippleUI.form(action = "/analysis/$uniqueID", method = "POST", [
      Stipple.select(:select, name="selectedResponse", options=:option,)

      btn("Analyze", type="submit", color="primary", @click(:process))
    ])

    card(class = "q-my-md", [
      card_section(h2("Analysis Information"))
      card_section(["Unique ID: ", uniqueID])
      card_section("Response Type: {{ output1 }}")
      card_section("{{ output2 }}")
    ])

  ]))
end

route("/") do
  model = @init()
  uniqueID = model.inputID[1:10]
  page(model, ui(uniqueID), title="BiDRA") |> html
end

route("/upload/:valID", method = POST) do
  uniqueID = Genie.Requests.payload(:valID)
  println(uniqueID)

  files = Genie.Requests.filespayload()
  for f in files
    write(joinpath(UPLOAD_PATH, "$uniqueID.csv"), f[2].data)
    @info "Uploading: " * f[2].name
  end

  if length(files) == 0
    @info "No file uploaded"
  end

end

route("/analysis/:valID", method=POST) do
  uniqueID = String(Genie.Requests.payload(:valID))
  select = Genie.Requests.payload(:selectedResponse)

  user_message = "Your experiments are currently being analyzed. Please do not reload nor close this page.\nOnce the process is completed, results will automatically be downloaded."
  user_message
  sleep(10)
  
  ### BiDRA Analysis process
  ## Define ResponseType + reload dataset
  user_message *= "\n\n-> Retrieving dataset and information: "
  user_message
  
  ANALYSIS_PATH = create_storage_dir("$RESULTS_PATH/$uniqueID")  
  responseType = getResponseTypeInt(select)
  dataset, expId, N = getDataset(UPLOAD_PATH, uniqueID)
  user_message *= " there are $N experiments to analyse."
  user_message

  ## Do inference
  user_message *= "\n\n-> Start inference for each experiment: [ "
  posterior_path = create_storage_dir("$ANALYSIS_PATH/posterior/")
  doseResponse_path = create_storage_dir("$ANALYSIS_PATH/doseResponse/")

  SD_ALL = []
  ΔWAICₖ_ALL = []
  posterior_ALL = DataFrame()

  for i in 1:N
      subsetId = expId[i]
    
      subsetDataset = filter(:id => x -> x == subsetId, dataset)
      posterior = doInference(subsetDataset, subsetId, responseType, posterior_path);
      posterior[!, :exp_id] = repeat([subsetId], size(posterior)[1])
      posterior_ALL =  vcat(posterior_ALL, posterior)
        
      plotDoseResponse(subsetDataset, posterior, subsetId, responseType, doseResponse_path);

      waicₖ_bidra, waicₖ_line = get_waicₖ(subsetDataset, posterior, responseType)
      ΔWAICₖ = waicₖ_bidra - waicₖ_line
      ΔWAICₖ_ALL = vcat(ΔWAICₖ_ALL, ΔWAICₖ)
      SD_ALL = vcat(SD_ALL, std(subsetDataset.response))
      user_message *= ". "
      user_message 
  end
  user_message *= "]"
  user_message

  ## Post-inference analysis
  dataset_metrics = DataFrame(exp_id=expId, SD=SD_ALL, ΔWAICₖ=ΔWAICₖ_ALL)
  dataset_metrics[!, :completeness] = map(r -> r ≥ 20 ? 1 : 2, dataset_metrics.SD)
  dataset_metrics[!, :WAICₖ_model] = map(r -> r < 0 ? 1 : 2, dataset_metrics.ΔWAICₖ)
  dataset_metrics[!, :group] = dataset_metrics.WAICₖ_model .+ dataset_metrics.completeness .- 1

  # Posterior informative potentiel
  plotGroupAssignation(dataset_metrics, ANALYSIS_PATH)
  user_message *= "\n\n-> Post-inference analysis: posterior informative potentiel  "
  user_message

  if N ≥ 3
    ## ranking
    rankingAnalysis(posterior_ALL, N, expId, 0, responseType, ANALYSIS_PATH)
    user_message *= " + IC50 rankings  "
    user_message

    rankingAnalysis(posterior_ALL, N, expId, 1, responseType, ANALYSIS_PATH)
    user_message *= " + HDR rankings  "
    user_message
  
  elseif  N == 2
    ## Comparison of two experiments
    plotPairedComparison(posterior_ALL, expId, ANALYSIS_PATH)
    user_message *= " + Comparison of two experiments  "
    user_message
  end

  user_message *= "\n\nANALYSIS COMPLETED -- Results will be automatically downloaded"
  user_message
  getPwd = pwd()
  cd(RESULTS_PATH)

  fn = createZip(uniqueID, getPwd, UPLOAD_PATH, ANALYSIS_PATH)
  HTTP.Response(200, ["Content-Type"=>"application/zip"], body=read("$RESULTS_PATH/$fn"))
  #user_message *= "\n\nYou will be redirected to the main page."
  
  #Genie.Renderer.redirect("/")
end

Genie.isrunning(:webserver) || up()
