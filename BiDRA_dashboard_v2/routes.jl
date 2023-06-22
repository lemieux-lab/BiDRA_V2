using Genie
using GenieSession
using Stipple, Stipple.ReactiveTools
using StippleUI
using HTTP

include("util.jl")
include("juPlot.jl")

@appname BiDRA_V2

const UPLOAD_PATH = create_storage_dir("Backend_Upload")
const RESULTS_PATH = create_storage_dir("Analysis_Results")

#GenieSession.__init__()

@app begin
  #println(:uniqueID)
  @in process = false
  @out inputID = getUniqueID()


  @out load = false
  @out output1 = "" ## response type selected by user
  @out output2 = "" ## analysis step 1 (retrieved data and info)
  @out output3 = "" ## analysis step 2 (do inference)
  @out output4 = "" ## analysis step 3 (Post-inference analysis: posterior informative potentiel)
  @out output5 = "" ## analysis completed
  @out output6 = "" ## new analysis info

  @in select = "Ascending"
  @in option = ["Ascending", "Descending"]

  @onbutton process begin
    println(input)
    load = true

    output1, output2, output3, output4, output5, output6 = resetOuput()
  
    output1 = select 

    ### BiDRA Analysis process
    ## Define ResponseType + reload dataset
    output2 = output2 * " Retrieved dataset and information: "
    ANALYSIS_PATH = create_storage_dir("$RESULTS_PATH/$inputID")
    
    responseType = getResponseTypeInt(select)
    dataset, expId, N = getDataset(UPLOAD_PATH, inputID)
    output2 = output2 * " there are $N experiments to analyse."

    ## Do inference
    output3 = output3 * "Start inference for each experiment: [ "
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
        output3 = output3 * ". "
    end
    output3 = output3 * "]"

    ## Post-inference analysis
    dataset_metrics = DataFrame(exp_id=expId, SD=SD_ALL, ΔWAICₖ=ΔWAICₖ_ALL)
    dataset_metrics[!, :completeness] = map(r -> r ≥ 20 ? 1 : 2, dataset_metrics.SD)
    dataset_metrics[!, :WAICₖ_model] = map(r -> r < 0 ? 1 : 2, dataset_metrics.ΔWAICₖ)
    dataset_metrics[!, :group] = dataset_metrics.WAICₖ_model .+ dataset_metrics.completeness .- 1

    # Posterior informative potentiel
    plotGroupAssignation(dataset_metrics, ANALYSIS_PATH)
    output4 = output4 * "Post-inference analysis: posterior informative potentiel  "

    if N ≥ 3
      ## ranking
      rankingAnalysis(posterior_ALL, N, expId, 0, responseType, ANALYSIS_PATH)
      output4 = output4 * " + IC50 rankings  "
      rankingAnalysis(posterior_ALL, N, expId, 1, responseType, ANALYSIS_PATH)
      output4 = output4 * " + HDR rankings  "
    elseif  N == 2
      ## Comparison of two experiments
      plotPairedComparison(posterior_ALL, expId, ANALYSIS_PATH)
      output4 = output4 * " + Comparison of two experiments  "
    end

    output5 = output5 * "ANALYSIS COMPLETED -- Results will be automatically downloaded"
    getPwd = pwd()
    cd(RESULTS_PATH)

    fn = createZip(inputID, getPwd, UPLOAD_PATH, RESULTS_PATH)
    HTTP.Response(200, ["Content-Type" => "application/zip"], body = read("$RESULTS_PATH/$fn"))
    output6 = output6 * "To do a new analysis, please reload the page."
  end
end

function ui(uniqueID)
  row(cell(class = "st-module", [
    uploader(name="fileUpload", label="Upload Dataset", accept=".csv", multiple=false, method="POST", url="/upload/$uniqueID", autoupload=true)

    Stipple.select(:select, name="selectedResponse", options=:option,)

    btn(class="q-my-md", name="process", "Analyze", @click(:process), color="green", disable=:load)

    card(class = "q-my-md", [
      card_section(h2("Analysis Information"))
      card_section("Response Type: {{ output1 }}")
      card_section(["Unique ID: ", uniqueID])
      card_section(["STARTING ANALYSIS PROCESS<br>{{ output2 }}<br>{{ output3 }}<br>{{ output4 }}<br><br>{{ output5 }}"])
      card_section(["{{ output6 }}"])
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

Genie.isrunning(:webserver) || up()
