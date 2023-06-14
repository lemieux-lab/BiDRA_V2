using Stipple, Stipple.ReactiveTools
using StippleUI

include("util.jl")

@appname Inverter

const FILE_PATH = create_storage_dir("Backend_Upload")

@app begin
  @in process = false
  @out output1 = ""
  @out output2 = ""
  @out output3 = ""

  @in select = "Ascending"
  @in option = ["Ascending", "Descending"]

  @onbutton process begin
    output1 = select
    output2 = uniqueID
    output3 = "Starting analysis process"
    
    responseType = getResponseTypeInt(select)

    output3 = output3 * "\nRetrieved dataset ... "
    sleep(3)
    output3 = output3 * "DONE (time)"

    sleep(3)
    N = 5
    output3 = output3 * "\n   --- There are $N individuel experiments to analyze "
    
    output3 = output3 * "\nStarting inference"
    for i in 1:N
        output3 = output3 * "\n   --- Exp. $(i): "
        sleep(3)
        output3 = output3 * "posterior inferred, "
        sleep(3)
        output3 = output3 * "DR curve plotted, "
        sleep(3)
        output3 = output3 * "informative potential calculated "
    end


  end
end

function ui()
    row(cell(class = "st-module", [
        uploader(label="Upload Dataset", accept=".csv", multiple=false, method="POST", url="/upload", field__name="csv_file")

        Stipple.select(:select, name="selectedResponse", options=:option,)

        btn(class = "q-my-md", "Analyze", @click(:process), color = "primary", loading="!process")
        
        card(class = "q-my-md", [
            card_section(h2("Analysis Information"))
            card_section("Response Type: {{ output1 }}")
            card_section(["Unique ID: {{ output2 }}"])
            card_section(["{{ output3 }}"])
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

        write(joinpath(FILE_PATH, "$uniqueID.csv"), f[2].data)
        @info "Uploading: " * f[2].name
    end
    if length(files) == 0
        @info "No file uploaded"
    end
    return "upload done"
end

Genie.isrunning(:webserver) || up()