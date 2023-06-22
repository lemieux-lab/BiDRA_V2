function llogistic(param::Array, responseType::Int) 
    HDR, LDR, ic50, slope = param 

    if responseType == 1
        return x -> LDR + ((HDR - LDR) / (1 + 10^(slope * (ic50 - x))))
    else
        return x -> HDR + ((LDR - HDR) / (1 + 10^(slope * (x - ic50))))
    end
end

@model function model_BIDRA(xs::Array, ys::Array, HDRmixture::MixtureModel, LDRμ::Int, responseType::Int) 
    HDR ~ HDRmixture
    LDR ~ Normal(LDRμ, 10)
    ic50 ~ Normal(0,10)
    slope ~ LogNormal(0.5, 1)
    
    σ ~ LogNormal(1, 1)

    for i in 1:length(xs)
        f = llogistic([HDR, LDR, ic50, slope], responseType)
        ys[i] ~ Normal(f(xs[i]), σ)
    end
end

@model function model_line(xs::Array, ys::Array, LDRμ::Int)
    line ~ Normal(LDRμ, 10)
    σ ~ LogNormal(1, 1)

    for i in 1:length(xs)
        ys[i] ~ Normal(line, σ)
    end
end