using Statistics, StatsBase
using CairoMakie

function plotDoseResponse(subsetData::DataFrame, posterior::DataFrame, id::String, responseType::Int, doseResponse_path::String)  
    concentration = subsetData.concentration
    response = subsetData.response 

    concentrationRange = collect(minimum(concentration)-1:0.1:maximum(concentration)+1)
    paramPosterior = posterior[:, [:HDR, :LDR, :ic50, :slope]]
    fitIte = map(row -> llogistic(Array(row), responseType), eachrow(paramPosterior))
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
    axislegend(axmain, position=responseType == 0 ? :lb : :lt,)

    ### IC50 posterior
    lineMed = round(median(posterior.ic50), digits=2)
    lineCI = round.(percentile(posterior.ic50, [2.5, 97.5]), digits=2)

    vspan!(axic50, lineCI[1], lineCI[2], color=(:orange, 0.3), label="95% C.I. $lineCI")
    hist!(axic50, posterior.ic50, bins=50, color=:black)
    vlines!(axic50, lineMed, color=:orange, label="Median $lineMed")

    if lineCI[2] > maximum(concentration)
        vlines!(axic50, [minimum(concentration), maximum(concentration)], color=:black, label="Exp. dose range")
    end
    
    axislegend(axic50, position = :lt,)
    hideydecorations!(axic50)

    ### LDR posterior
    lineMed = round(median(posterior.LDR), digits=2)
    lineCI = round.(percentile(posterior.LDR, [2.5, 97.5]), digits=2)
    
    hspan!(axldr, lineCI[1], lineCI[2], color=(:blue, 0.3), label="95% C.I. LDR\n$lineCI")
    hist!(axldr, posterior.LDR, bins=50, color=:black, direction=:x)
    hlines!(axldr, lineMed, color=:blue, label="Median LDR\n$lineMed")

    linkyaxes!(axmain, axldr)
    hideydecorations!(axldr, grid=false)
    hidexdecorations!(axldr, label=false)
    Legend(f[1, 4], axldr)

    ### HDR posterior
    lineMed = round(median(posterior.HDR), digits=2)
    lineCI = round.(percentile(posterior.HDR, [2.5, 97.5]), digits=2)

    hspan!(axhdr, lineCI[1], lineCI[2], color=(:green, 0.3), label="95% C.I. HDR\n$lineCI")
    hist!(axhdr, posterior.HDR, bins=50, color=:black, direction=:x)
    hlines!(axhdr, lineMed, color=:green, label="Median HDR\n$lineMed")

    linkyaxes!(axmain, axhdr)
    hideydecorations!(axhdr, grid=false)
    hidexdecorations!(axhdr, label=false)
    Legend(f[2, 4], axhdr)

    ### Slope posterior
    lineMed = round(median(posterior.slope), digits=2)
    lineCI = round.(percentile(posterior.slope, [2.5, 97.5]), digits=2)

    vspan!(axslope, lineCI[1], lineCI[2], color=(:purple, 0.3), label="95% C.I. Slope\n$lineCI")
    hist!(axslope, posterior.slope, bins=50, color=:black)
    vlines!(axslope, lineMed, color=:purple, label="Median Slope\n$lineMed")
    
    Legend(f[3, 4], axslope)
    hideydecorations!(axslope)

    ### σ posterior
    lineMed = round(median(posterior.σ), digits=2)
    lineCI = round.(percentile(posterior.σ, [2.5, 97.5]), digits=2)

    vspan!(axσ, lineCI[1], lineCI[2], color=(:pink, 0.3), label="95% C.I. σ\n$lineCI")
    hist!(axσ, posterior.σ, bins=50, color=:black)
    vlines!(axσ, lineMed, color=:pink, label="Median σ\n$lineMed")

    Legend(f[4, 4], axσ)
    hideydecorations!(axσ)

    ### Save Figure
    fn = "doseResponse_expID_$id"*".pdf"
    CairoMakie.save("$doseResponse_path$fn", f)
end

function plotGroupAssignation(dataset_metrics::DataFrame, ANALYSIS_PATH::String)
    N = unique(dataset_metrics.exp_id) |> length
    w = 250*3
    h = 250+(25*N)

    colorVal = ["green", "orange", "red"]

    f = Figure(backgroundcolor="transparent", resolution=(w, h))
    ax = Axis(f[1, 1], title="Posterior Informative Potentiel Flag")

    scatter!(ax, repeat([1], N), collect(1:N), color=colorVal[dataset_metrics.group], markersize=10)
    text!(ax, repeat([1.5], N), collect(1:N), text=dataset_metrics.exp_id, align=(:left, :center))
    
    xlims!(ax, 0.5, 10)
    hidedecorations!(ax)
    hidespines!(ax)

    fn = "informative_potential_posterior"*".pdf"
    CairoMakie.save("$ANALYSIS_PATH/$fn", f)
end

function plotPairedComparison(posterior::DataFrame, expId::Vector{String}, ANALYSIS_PATH::String)
    posterior_1 = filter(:exp_id => x -> x == expId[1], posterior)
    posterior_2 = filter(:exp_id => x -> x == expId[2], posterior)
    
    function plotCol(c::Int, p::Symbol)
        ax_up = Axis(f[1, c], title=String(p), xlabel="Δ posterior", ylabel="Count")
        ax_low = Axis(f[2, c], title="$(String(p)) QQ posterior", xlabel="Exp. $(expId[1])", ylabel="Exp. $(expId[2])")

        Δ = posterior_1[:, p] .- posterior_2[:, p]
        pct_1 = round(length(Δ[Δ .>= 0]) / length(Δ), digits=2)
        pct_2 = round(length(Δ[Δ .< 0]) / length(Δ), digits=2)

        hist!(ax_up, Δ, bins=50, color=:black)
        vlines!(ax_up, [0.], color=:blue, linestyle=:dash)
        text!(ax_up, 1., 0.5, text="Exp. $(expId[1])\nis larger ($pct_1)", space=:relative, rotation=3*pi/2, align=[:center, :top], color=:blue)
        text!(ax_up, 0., 0.5, text="Exp. $(expId[2])\nis larger ($pct_2)", space=:relative, rotation=pi/2, align=[:center, :top], color=:blue)

        scatter!(ax_low, sort(posterior_1[:, p]), sort(posterior_2[:, p]), color=(:black, 0.1))
        ablines!(ax_low, [0], [1], color=:blue, linestyle=:dash)

        return ax_up, ax_low
    end

    f = Figure(backgroundcolor="transparent", resolution=(300*5, 600))

    axA_1, axA_2 = plotCol(1, :ic50)
    axB_1, axB_2 = plotCol(2, :HDR)
    axC_1, axC_2 = plotCol(3, :LDR)
    axD_1, axD_2 = plotCol(4, :slope)
    axE_1, axE_2 = plotCol(5, :σ)

    fn = "posterior_comparison"*".pdf"
    CairoMakie.save("$ANALYSIS_PATH/$fn", f)
end

function plotRanking(sorted_rank_matrix::Matrix{Float64}, posterior::DataFrame, sorted_expId::Vector, expId::Vector{String}, N::Int, n_samples::Int, paramRank::Int, ANALYSIS_PATH::String)
    function getPosteriorBinRange(posterior::DataFrame, pr::Symbol, pr_bins::Vector, expId_order::Vector, N::Int, n_samples::Int)
        ## Matrix of posterior samplings (sampling x compound)
        posterior_mtx = reshape(posterior[:, pr], (n_samples, N))
      
        ## Calcute histogram binning (bins x compound)
        hist_mtx = mapreduce(c -> fit(Histogram, c, pr_bins).weights, hcat, eachcol(posterior_mtx))
      
        ## Sort by ranking order 
        hist_sorted = hist_mtx[:, expId_order]
        return hist_sorted
    end
      
    w = N*10 + (400*2)
    h = N*30 + 60

    ## info for posterior representations
    eff_metrics = [:ic50, :HDR]
    minIC50, maxIC50 = percentile(posterior.ic50, [2.5, 97.5])#extrema(posterior.ic50)
    minHDR, maxHDR = percentile(posterior.HDR, [2.5, 97.5])#extrema(posterior.HDR)
    bins_range = [collect(floor(minIC50):0.1:ceil(maxIC50)), collect(floor(minHDR):1:ceil(maxHDR))]

    binRange_ic50 = getPosteriorBinRange(posterior, eff_metrics[1], bins_range[1], sorted_expId, N, n_samples)
    binRange_hdr = getPosteriorBinRange(posterior, eff_metrics[2], bins_range[2], sorted_expId, N, n_samples)

    xticksIC50 = range(1, length(bins_range[1]), 4)
    xticksHDR = range(1, length(bins_range[2]), 4)

    if paramRank == 0
        PARAM = "IC50"
    elseif paramRank == 1
        PARAM = "HDR"
    end
      
    f = Figure(backgroundcolor="transparent", resolution=(w, h));
    ax = Axis(f[1, 1], xlabel="$PARAM Ranks", xticks=(1:5:N, string.(1:5:N)), yticks=(1:N, expId[sorted_expId]))
      
    hm = CairoMakie.heatmap!(ax, sorted_rank_matrix, colormap=["gray95", "gray35", "black"])
    Colorbar(f[2, 1], hm, label = "Rank Prob.", vertical = false, flipaxis = false)
      
    axIc50 = Axis(f[1, 2], xlabel="IC50 posterior", xticks=(xticksIC50, string.(bins_range[1][Int.(round.(collect(xticksIC50)))])))
    hm = CairoMakie.heatmap!(axIc50, binRange_ic50, colormap=["gray95", "gray35", "black"])
    hideydecorations!(axIc50)
    Colorbar(f[2, 2], hm, label="Rank Prob.", vertical = false, flipaxis = false)
      
    axHDR = Axis(f[1, 3], xlabel="HDR posterior", xticks=(xticksHDR, string.(bins_range[2][Int.(round.(collect(xticksHDR)))])))
    hm = CairoMakie.heatmap!(axHDR, binRange_hdr, colormap=["gray95", "gray35", "black"])
    Colorbar(f[2, 3], hm, label="Count", vertical = false, flipaxis = false)
    hideydecorations!(axHDR)
      
    colsize!(f.layout, 1, Relative(2/3))
    linkyaxes!(ax, axIc50, axHDR)

    fn = "$(PARAM)_rankings"*".pdf"
    CairoMakie.save("$ANALYSIS_PATH/$fn", f)
end

function plotPairWiseParameters(posterior::DataFrame)
    f = Figure(backgroundcolor="transparent", resolution=(300*3, 300*3))
    axA = Axis(f[1, 1], xlabel="IC₅₀", ylabel="HDR")
    axB = Axis(f[2, 1], xlabel="IC₅₀", ylabel="Slope")
    axC = Axis(f[3, 1], xlabel="IC₅₀", ylabel="LDR")
    axD = Axis(f[2, 2], xlabel="HDR", ylabel="Slope")
    axE = Axis(f[3, 2], xlabel="HDR", ylabel="LDR")
    axF = Axis(f[3, 3], xlabel="Slope", ylabel="LDR")

    scatter!(axA, posterior.ic50, posterior.HDR, color=(:blue, 0.1))
    scatter!(axA, sort(posterior.ic50), sort(posterior.HDR), color=(:black, 0.1))
    scatter!(axB, posterior.ic50, posterior.slope, color=(:blue, 0.1))
    scatter!(axB, sort(posterior.ic50), sort(posterior.slope), color=(:black, 0.1))
    scatter!(axC, posterior.ic50, posterior.LDR, color=(:blue, 0.1))
    scatter!(axC, sort(posterior.ic50), sort(posterior.LDR), color=(:black, 0.1))
    scatter!(axD, posterior.HDR, posterior.slope, color=(:blue, 0.1))
    scatter!(axD, sort(posterior.HDR), sort(posterior.slope), color=(:black, 0.1))
    scatter!(axE, posterior.HDR, posterior.LDR, color=(:blue, 0.1))
    scatter!(axE, sort(posterior.HDR), sort(posterior.LDR), color=(:black, 0.1))
    scatter!(axF, posterior.slope, posterior.LDR, color=(:blue, 0.1))
    scatter!(axF, sort(posterior.slope), sort(posterior.LDR), color=(:black, 0.1))

    linkxaxes!(axA, axB, axC)
    linkxaxes!(axD, axE)
    linkyaxes!(axB, axD)
    linkyaxes!(axC, axE, axF)

    hidexdecorations!(axA, grid=false)
    hidexdecorations!(axB, grid=false)
    hidedecorations!(axD, greid=false)
    hideydecorations!(axE, grid=false)
    hideydecorations!(axF, grid=false)
end
