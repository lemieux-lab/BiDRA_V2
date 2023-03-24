using CairoMakie

function plotDoseResponse(subsetData::DataFrame, posterior::DataFrame, id::String)  
    concentration = subsetData.concentration
    response = subsetData.response 

    concentrationRange = collect(minimum(concentration)-1:0.1:maximum(concentration)+1)
    paramPosterior = posterior[:, [:HDR, :LDR, :ic50, :slope]]
    fitIte = map(row -> llogistic(Array(row)), eachrow(paramPosterior))
    responsePosterior = mapreduce(f -> f.(concentrationRange), hcat, fitIte)

    medianCurve = median.(eachrow(responsePosterior))
    upperBound = percentile.(eachrow(responsePosterior), 97.5)
    lowerBound = percentile.(eachrow(responsePosterior), 2.5)

    f = Figure(resolution=(1200, 800), backgroundcolor=:transparent)
    axmain = Axis(f[1:2, 1], width=600, title="ExpId $id", xlabel="log10 Concentrations", ylabel="Response (%)")
    axic50 = Axis(f[3:4, 1], xlabel="IC₅₀ Posterior")
    axldr = Axis(f[1:2, 2], width=120, xlabel="LDR Posterior")
    axhdr = Axis(f[1:2, 3], width=120, xlabel="HDR Posterior")
    axslope = Axis(f[3, 2:3], height=120, xlabel="Slope Posterior")
    axσ = Axis(f[4, 2:3],height=120, xlabel="σ Posterior")

    ### Main plot - dose response curve
    band!(axmain, concentrationRange, lowerBound, upperBound, color=(:gray, 0.4), label="95% C.I.")
    scatter!(axmain, concentration, response, color=:black, label="Exp. Responses")
    lines!(axmain, concentrationRange, medianCurve, color=:black, label="Median curve")
    axislegend(axmain, position = :lt,)

    ### IC50 posterior
    hist!(axic50, posterior.ic50, bins=50, color=:black)
    if maximum(posterior.ic50) > maximum(concentration)
        vspan!(axic50, minimum(concentration), maximum(concentration), color=(:orange, 0.2), label="Exp. dose range")
    end
    lineMed = round(median(posterior.ic50), digits=2)
    lineCI = round.(percentile(posterior.ic50, [2.5, 97.5]), digits=2)
    vlines!(axic50, lineMed, color=:orange, label="Median $lineMed")
    vlines!(axic50, lineCI, color=:orange, linestyle=:dash, label="95% C.I. $lineCI")
    axislegend(axic50, position = :rt,)
    hideydecorations!(axic50)

    ### LDR posterior
    lineMed = round(median(posterior.LDR), digits=2)
    lineCI = round.(percentile(posterior.LDR, [2.5, 97.5]), digits=2)
    hist!(axldr, posterior.LDR, bins=50, color=:black, direction=:x)
    hlines!(axldr, lineMed, color=:blue, label="Median LDR\n$lineMed")
    hlines!(axldr, lineCI, color=:blue, linestyle=:dash, label="95% C.I. LDR\n$lineCI")
    linkyaxes!(axmain, axldr)
    hideydecorations!(axldr, grid=false)
    hidexdecorations!(axldr, label=false)
    Legend(f[1, 4], axldr)

    ### HDR posterior
    lineMed = round(median(posterior.HDR), digits=2)
    lineCI = round.(percentile(posterior.HDR, [2.5, 97.5]), digits=2)
    hist!(axhdr, posterior.HDR, bins=50, color=:black, direction=:x)
    hlines!(axhdr, lineMed, color=:green, label="Median HDR\n$lineMed")
    hlines!(axhdr, lineCI, color=:green, linestyle=:dash, label="95% C.I. HDR\n$lineCI")
    linkyaxes!(axmain, axhdr)
    hideydecorations!(axhdr, grid=false)
    hidexdecorations!(axhdr, label=false)
    Legend(f[2, 4], axhdr)

    ### Slope posterior
    lineMed = round(median(posterior.slope), digits=2)
    lineCI = round.(percentile(posterior.slope, [2.5, 97.5]), digits=2)
    hist!(axslope, posterior.slope, bins=50, color=:black)
    vlines!(axslope, lineMed, color=:purple, label="Median Slope\n$lineMed")
    vlines!(axslope, lineCI, color=:purple, linestyle=:dash, label="95% C.I. Slope\n$lineCI")
    Legend(f[3, 4], axslope)
    hideydecorations!(axslope)

    ###
    lineMed = round(median(posterior.σₛ), digits=2)
    lineCI = round.(percentile(posterior.σₛ, [2.5, 97.5]), digits=2)
    hist!(axσ, posterior.σₛ, bins=50, color=:black)
    vlines!(axσ, lineMed, color=:pink, label="Median σ\n$lineMed")
    vlines!(axσ, lineCI, color=:pink, linestyle=:dash, label="95% C.I. σ\n$lineCI")
    Legend(f[4, 4], axσ)
    hideydecorations!(axσ)

    doseResponse_path = dataset_path*"doseResponse/"
    mkpath(doseResponse_path)
    
    fn = "doseResponse_expID_$id"*".pdf"
    CairoMakie.save("$doseResponse_path$fn", f)
    #fn = "doseResponse_expID_$id"*".png"
    #CairoMakie.save("$mainPath/figure/$fn", f)
end