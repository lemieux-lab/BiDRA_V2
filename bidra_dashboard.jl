using Stipple, Stipple.ReactiveTools
using StippleUI
using HTTP

include("util.jl")
include("juPlot.jl")

@appname BiDRA_V2

const UPLOAD_PATH = create_storage_dir("Backend_Upload")
const RESULTS_PATH = create_storage_dir("Analysis_Results")

@app begin
  @in process = false

  @out load = false
  @out output1 = "" ## response type selected by user
  @out output2 = "" ## uniqueID
  @out output3 = "" ## analysis step 1 (retrieved data and info)
  @out output4 = "" ## analysis step 2 (do inference)
  @out output5 = "" ## analysis step 3 (Post-inference analysis: posterior informative potentiel)
  @out output6 = "" ## analysis completed

  @in select = "Ascending"
  @in option = ["Ascending", "Descending"]

  @onbutton process begin
    load = true

    output1, output2, output3, output4, output5, output6 = resetOuput()
  
    output1 = select 
    output2 = uniqueID 

    ### BiDRA Analysis process
    ## Define ResponseType + reload dataset
    output3 = output3 * " Retrieved dataset and information: "
    global ANALYSIS_PATH = create_storage_dir("$RESULTS_PATH/$uniqueID")
    responseType = getResponseTypeInt(select)
    dataset, expId, N = getDataset(UPLOAD_PATH, uniqueID)
    output3 = output3 * " there are $N experiments to analyse."

    ## Do inference
    output4 = output4 * "Start inference for each experiment: [ "
    global posterior_path = create_storage_dir("$ANALYSIS_PATH/posterior/")
    global doseResponse_path = create_storage_dir("$ANALYSIS_PATH/doseResponse/")

    SD_ALL = []
    ΔWAICₖ_ALL = []
    posterior_ALL = DataFrame()

    for i in 1:N
        subsetId = expId[i]
    
        subsetDataset = filter(:id => x -> x == subsetId, dataset)
        posterior = doInference(subsetDataset, subsetId, responseType);
        posterior[!, :exp_id] = repeat([subsetId], size(posterior)[1])
        posterior_ALL =  vcat(posterior_ALL, posterior)
        
        plotDoseResponse(subsetDataset, posterior, subsetId, responseType);

        waicₖ_bidra, waicₖ_line = get_waicₖ(subsetDataset, posterior, responseType)
        ΔWAICₖ = waicₖ_bidra - waicₖ_line
        ΔWAICₖ_ALL = vcat(ΔWAICₖ_ALL, ΔWAICₖ)
        SD_ALL = vcat(SD_ALL, std(subsetDataset.response))
        output4 = output4 * ". "
    end
    output4 = output4 * "]"

    ## Post-inference analysis
    dataset_metrics = DataFrame(exp_id=expId, SD=SD_ALL, ΔWAICₖ=ΔWAICₖ_ALL)
    dataset_metrics[!, :completeness] = map(r -> r ≥ 20 ? 1 : 2, dataset_metrics.SD)
    dataset_metrics[!, :WAICₖ_model] = map(r -> r < 0 ? 1 : 2, dataset_metrics.ΔWAICₖ)
    dataset_metrics[!, :group] = dataset_metrics.WAICₖ_model .+ dataset_metrics.completeness .- 1

    # Posterior informative potentiel
    plotGroupAssignation(dataset_metrics)
    output5 = output5 * "Post-inference analysis: posterior informative potentiel  "

    if N ≥ 3
      ## ranking
      rankingAnalysis(posterior_ALL, N, expId, 0, responseType)
      output5 = output5 * " + IC50 rankings  "
      rankingAnalysis(posterior_ALL, N, expId, 1, responseType)
      output5 = output5 * " + HDR rankings  "
    elseif  N == 2
      ## Comparison of two experiments
      plotPairedComparison(posterior_ALL, expId)
      output5 = output5 * " + Comparison of two experiments  "
    end

    output6 = output6 * "ANALYSIS COMPLETED -- Results will be automatically downloaded"
    getPwd = pwd()
    cd(RESULTS_PATH)

    fn = createZip(uniqueID, getPwd)
    HTTP.Response(200, ["Content-Type" => "application/zip"], body = read("$RESULTS_PATH/$fn"))

    global uniqueID = ""
    load = false
  end
end

function ui()
  row(cell(class = "st-module", [
    uploader(label="Upload Dataset", accept=".csv", multiple=false, method="POST", url="/upload", field__name="csv_file", autoupload=true)

    Stipple.select(:select, name="selectedResponse", options=:option,)

    btn(class = "q-my-md", name="process", "Analyze", @click(:process), color = "green", disable=:load)

    card(class = "q-my-md", [
      card_section(h2("Analysis Information"))
      card_section("Response Type: {{ output1 }}")
      card_section(["Unique ID: {{ output2 }}"])
      card_section(["STARTING ANALYSIS PROCESS<br>{{ output3 }}<br>{{ output4 }}<br>{{ output5 }}<br><br>{{ output6 }}"])
      #card_section([btn(class = "q-my-md", name="new", "New Analysis", @click(:new), color = "primary", )])
    ])
  ]))
end

route("/") do
  model = @init
  page(model, ui(), title="BiDRA") |> html
end

route("/upload", method = POST) do
    files = Genie.Requests.filespayload()
    for f in files
        global uniqueID = getUniqueID()
        global experimentName = split(f[2].name, ".")[1]

        write(joinpath(UPLOAD_PATH, "$uniqueID.csv"), f[2].data)
        @info "Uploading: " * f[2].name
    end
    if length(files) == 0
        @info "No file uploaded"
    end
    return "upload done"
end

Genie.isrunning(:webserver) || up()
