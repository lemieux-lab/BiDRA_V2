using Genie, Stipple
using Genie.Requests
using StippleUI

function create_storage_dir(name)
    try
      mkdir(joinpath(@__DIR__, name)) 
    catch 
      @warn "directory already exists" 
    end
    return joinpath(@__DIR__, name)
end

# Generate file path
const FILE_PATH = create_storage_dir("Backend_Upload")

# Define react model
@vars APP begin
    process::R{Bool} = false
    responseSelected::R{String} = "Ascending"
    responseOptions::R{Vector{String}} = ["Ascending", "Descending"]
end

function ui(model::APP)
    page(model, title="Dashboard",
    [
        heading("Dashboard") 

        row([
          Html.div(class="col-md-12", [
            uploader(label="Upload Dataset", accept=".csv", multiple=false, method="POST", url="/upload", field__name="csv_file")
          ])
        ])

        row([
          Html.div(class="col-md-4", [
            Stipple.select(:responseSelected, name="selectedResponse", options=:responseOptions,)
          ])
        ])

        row([
          Html.div(class="col-md-4", [
            btn(label="Analyze", color="secondary", @click("process=true"))
          ])
        ])
    ])
  end

route("/") do
    APP |> init |> ui |> html

    if process
      println("Start sleep")
      sleep(15)
      println("end sleep")
      process = false
    end
end
  
#uploading csv files to the backend server
route("/upload", method = POST) do
    files = Genie.Requests.filespayload()
    for f in files
        write(joinpath(FILE_PATH, f[2].name), f[2].data)
        @info "Uploading: " * f[2].name
    end
    if length(files) == 0
        @info "No file uploaded"
    end
    return "upload done"
end

up()
