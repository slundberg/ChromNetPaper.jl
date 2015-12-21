using MLBase
using JSON
using SimplePlot

export
    network_enrichment,
    id2uniprot,
    id2truth,
    uniprot2truth,
    id2celltype,
    id2treatments,
    id2target,
    ishistone,
    random_cor_matrix,
    conditional_cov,
    truth_matrix,
    mask_matrix,
    unique_ppi_pairs,
    apply_to_celltypes,
    network_enrichment_rank,
    bootstrap_network_enrichment_rank,
    network_enrichment_density,
    enrichment_rank,
    bin_values_avg

"Simple helper method to get the upper off-diagonal entries of a matrix."
function upper{T}(X::Array{T,2})
    x = T[]
    for i in 1:size(X)[2], j in 1:i-1
        push!(x, X[j,i])
    end
    x
end

"""
This reports the enrichment (over random) of true edges from the given score matrix.

X - The matrix of scores (higher is stronger)
T - A matrix of truth pairings.
M - An optional boolean matrix masking which entries should be evaluated.
"""
function network_enrichment(X::AbstractMatrix, T::Array{Bool,2}, M::Array{Bool,2}; weights=nothing, numEdges=:num_true)
    weights = weights == nothing ? ones(size(X)...) : weights
    @assert size(X) == size(T)
    @assert size(X) == size(M)
    @assert size(X) == size(weights)

    mask = upper(M)
    scores = upper(X)[mask]
    truth = upper(T)[mask]
    weightsMasked = upper(weights)[mask]
    if numEdges == :num_true
        numEdges = sum(weightsMasked[truth])
    end

    totalCorrect = 0.0
    totalWeight = 0.0
    sortedInds = sortperm(-scores)
    for ind in sortedInds
        totalWeight += weightsMasked[ind]
        if truth[ind]
            totalCorrect += weightsMasked[ind]
        end
        if totalWeight >= numEdges
            break
        end
    end
    (totalCorrect/totalWeight) / (sum(weightsMasked[truth])/sum(weightsMasked))
end
function network_enrichment(X::AbstractMatrix, T::Array{Bool,2}; weights=nothing, numEdges=:num_true)
    network_enrichment(X, T, ones(Bool, size(X)...); weights=weights, numEdges=numEdges)
end

function bootstrap_network_enrichment_rank(data, T, ids; numSamples=10, samplesAlpha=0.05, density=false, numBins=100)
    targets = unique([id2target(id) for id in filter(id->!ishistone(id), ids)])

    sendto(workers(), ChromNetPaper, sentData=data, sentT=T)
    samples = pmap(Progress(numSamples), i->begin

        # get a bootstrap re-sample
        bootstrapSample = StatsBase.sample(targets, length(targets))
        idWeights = Dict{ASCIIString,Float64}()
        for id in bootstrapSample
            idWeights[id] = get(idWeights, id, 0) + 1
        end

        # build a bootstrap weighting matrix
        w = zeros(length(ids))
        for j in 1:length(ids)
            w[j] = get(idWeights, ChromNetPaper.id2target(ids[j]), 0)
        end
        W = w*w';

        datum = Any[]
        for d in sentData
            x,y = ChromNetPaper.network_enrichment_rank(abs(d[1]), sentT, d[2], weights=W)
            push!(datum, (x,y))
        end
        datum
    end, 1:numSamples)
    sleep(0.1) # let the progress meter finish printing

    # bin the curves then line them up and average them
    numBins = numBins
    maxRank = minimum([round(Int, sum(d[2]) / 2) for d in data]) # find the smallest number of considered edges
    resSamples = Any[]
    resAvg = Any[]
    resAucs = Any[]
    for d in data
        push!(resSamples, zeros(numBins, numSamples))
        push!(resAvg, zeros(numBins))
        push!(resAucs, zeros(numSamples))
    end
    for i in 1:numSamples
        for j in 1:length(data)
            resSamples[j][:,i] = bin_values_avg(samples[i][j][1][1:maxRank], samples[i][j][2][1:maxRank], numBins)
            resAvg[j] .+= resSamples[j][:,i]
            resAucs[j][i] = sum(resSamples[j][:,i])
        end
    end
    for j in 1:length(data)
        resAvg[j] ./= numSamples
    end

    # scale the xs depending on if we are plotting by rank or density (of the smallest network given)
    optionalArgs = Dict()
    xs = linspace(1,maxRank,numBins)
    if density
        xs = linspace(0,1,numBins)
        optionalArgs[:xticks] = [0, 0.5, 1]
        optionalArgs[:xticklabels] = ["0%", "50%", "100%"]
    end

    layers = [
        [line(xs, resAvg[j], color=SimplePlot.defaultColors[j], data[j][3], linewidth=3) for j in 1:length(data)];
        vcat([[line(
                xs,
                resSamples[j][:,i],
                color=SimplePlot.defaultColors[j],
                alpha=samplesAlpha
            ) for i in 1:min(numSamples,100)] for j in 1:length(data)]...)
    ]

    layers,resAucs,maxRank
end

# fill in any blank spots with nearby values
function fill_nans(vals)
    newVals = copy(vals)
    for i in 2:length(vals)
        if isnan(vals[i])
            newVals[i] = newVals[i-1]
        end
    end
    for i in length(vals)-1:-1:1
        if isnan(vals[i])
            newVals[i] = newVals[i+1]
        end
    end
    newVals
end

# compute the rank enrichment curve
function enrichment_rank(truth, pred; weights=nothing)
    weights = weights == nothing ? ones(length(truth)) : weights
    p = sortperm(pred, rev=true)
    enrichment_rank(truth[p], weights=weights[p])
end
function enrichment_rank(rankedTruth; weights=nothing)
    weights = weights == nothing ? ones(length(rankedTruth)) : weights
    x = zeros(Float64, length(rankedTruth))
    y = zeros(Float64, length(rankedTruth))
    randRate = dot(rankedTruth,weights) / sum(weights)
    totalMatchedWeight = 0
    totalSeenWeight = 0

    for i in 1:length(rankedTruth)
        totalMatchedWeight += rankedTruth[i] * weights[i]
        totalSeenWeight += weights[i]
        x[i] = totalSeenWeight
        y[i] = (totalMatchedWeight/totalSeenWeight) / randRate
    end
    x,fill_nans(y)
end

function network_enrichment_rank(X::AbstractMatrix, T::Array{Bool,2}, M::Array{Bool,2}; weights=nothing)
    weights = weights == nothing ? ones(size(X)...) : weights
    @assert size(X) == size(T)
    @assert size(X) == size(M)
    @assert size(X) == size(weights)

    mask = upper(M)
    scores = upper(X)[mask]
    truth = upper(T)[mask]
    weightsMasked = upper(weights)[mask]
    x,y = enrichment_rank(truth, scores; weights=weightsMasked)
    x, y
end

function network_enrichment_density(X::AbstractMatrix, T::Array{Bool,2}, M::Array{Bool,2})
    x,y = network_enrichment_rank(X, T, M)
    x ./ x[end], y
end

uniprotHistones = ["Q71DI3", "P0C0S5", "P62805", "P84243"]
function ishistone(id)
    return id2uniprot(id) in uniprotHistones
end


"Load the matrix of pairwise dataset connections supported by BioGRID"
function truth_matrix(header::AbstractArray)
    P = length(header)
    T = zeros(Bool, P,P)
    for i in 1:P, j in i:P
        T[j,i] = T[i,j] = id2truth(header[j], header[i])
    end
    T
end

"Apply the passed function to each cell type sub-block of the given matrix."
function apply_to_celltypes(f, C, ids; cellTypes=nothing)
    out = zeros(size(C)...)
    if cellTypes == nothing
        cellTypes = unique([ChromNetPaper.id2celltype(id) for id in ids])
    end
    for cellType in cellTypes
        inds = find([ChromNetPaper.id2celltype(id) == cellType for id in ids])
        out[inds,inds] = f(C[inds,inds], ids[inds])
    end
    out
end

id2truthDict = Dict()
"Check if two datasets are connected in BioGRID."
function id2truth(id1, id2)
    uniprot2truth(id2uniprot(id1), id2uniprot(id2))
end
function uniprot2truth(uid1, uid2)
    global id2truthDict

    if (uid1 == uid2) && (uid1 != "")
        return true
    end

    # lazy load the data
    if length(id2truthDict) == 0
        f = open(joinpath(globalDataDir, "biogrid_human_swissprot.txt"))
        for line in eachline(f)
            parts = split(strip(line), '\t')
            id2truthDict[(parts[1],parts[2])] = true
            id2truthDict[(parts[2],parts[1])] = true
        end
        close(f)
    end

    get(id2truthDict, (uid1, uid2), false)
end

metadataDict = open(f->JSON.parse(readall(f)), joinpath(globalDataDir, "metadata.json"))

"Get the uniprot id of a given experiment id."
function id2uniprot(id)
    if haskey(metadataDict, id)
        return metadataDict[id]["targetUniProt"]
    end
    ""
end

"Get the cell type of a given experiment id."
function id2celltype(id)
    if haskey(metadataDict, id)
        return metadataDict[id]["cellType"]
    end
    ""
end

"Get the target name of a given experiment id."
function id2target(id)
    if haskey(metadataDict, id)
        return metadataDict[id]["name"]
    end
    ""
end

"Get the treatments of a given experiment id."
function id2treatments(id)
    if haskey(metadataDict, id)
        return metadataDict[id]["treatments"]
    end
    ""
end

# from: http://stackoverflow.com/questions/27677399/julia-how-to-copy-data-to-another-processor-in-julia
function sendto(p::Int, m=Main; args...)
    for (nm, val) in args
        @spawnat(p, eval(m, Expr(:(=), nm, val)))
    end
end
function sendto(ps::Vector{Int}, m=Main; args...)
    for p in ps
        sendto(p, m; args...)
    end
end

function mask_matrix(maskType, ids; excludeCellType=nothing, includeCrossEdges=false, histoneProteinLinks=false)
    P = length(ids)
    M = zeros(Bool, P, P)
    for i in 1:P, j in i:P

        # exclude connections involving histones
        exclude = false
        if histoneProteinLinks
            # only consider links involving one histone
            exclude = !ishistone(ids[i]) && !ishistone(ids[j])
            exclude = exclude || (ishistone(ids[i]) && ishistone(ids[j]))
        else
            # exclude connections involving histones
            exclude = ishistone(ids[i]) || ishistone(ids[j])
        end

        if excludeCellType != nothing
            exclude = exclude || id2celltype(ids[i]) == excludeCellType || id2celltype(ids[j]) == excludeCellType
        end

        # exclude same uniprot target connections
        exclude = exclude || (id2uniprot(ids[i]) == id2uniprot(ids[j]))

        if maskType == "within_all" || maskType == "within_all_separate"
            exclude = exclude || id2celltype(ids[i]) != id2celltype(ids[j])
        elseif maskType == "between_all"
            exclude = exclude || id2celltype(ids[i]) == id2celltype(ids[j])
        elseif maskType == "all"
            # no filter
        else # restrict by cellType
            if includeCrossEdges
                exclude = exclude || (id2celltype(ids[i]) != maskType && id2celltype(ids[j]) != maskType)
            else
                exclude = exclude || id2celltype(ids[i]) != maskType
                exclude = exclude || id2celltype(ids[j]) != maskType
            end
        end

        M[i,j] = M[j,i] = !exclude
    end
    M
end

function bin_values_avg(x, y, nbins)
    @assert length(x) == length(y)
    counts = zeros(nbins)
    sums = zeros(nbins)
    xmax = maximum(x)

    # bin the data
    for i in 1:length(x)
        ind = min(round(Int, floor(x[i]/xmax * nbins))+1, nbins)
        counts[ind] += 1
        sums[ind] += y[i]
    end

    ChromNetPaper.fill_nans(sums ./ counts)
end

function random_cor_matrix(K, netDensity)
    IC = diagm(abs(randn(K)))
    for i in 1:K, j in 1:i-1
        IC[i,j] = rand() < netDensity ? rand()-1.01 : 0.0
        IC[j,i] = IC[i,j]
    end
    mineval = minimum(eig(IC)[1])
    if mineval < 0
        IC -= eye(K)*mineval*1.01
    end
    C = inv(IC)
    Base.cov2cor!(C, sqrt(diag(C)))
end

function conditional_cov(C, controlInds, usedInds, reg=1e-8)
    controlC = C[controlInds,controlInds]
    targetC = C[usedInds,usedInds]
    crossC = C[usedInds,controlInds]
    A = targetC .- crossC*inv(controlC .+ reg*eye(size(controlC)[1]))*transpose(crossC) .+ reg*eye(size(targetC)[1])
end

function unique_ppi_pairs(X::AbstractMatrix, ids, T::Array{Bool,2}, M::Array{Bool,2}; numEdges=nothing, scoreThreshold=nothing)
    @assert size(X)[1] == length(ids)
    @assert size(X) == size(M)

    pairsMatrix = Array(Any, size(X)...)
    for i in 1:size(X)[1], j in i+1:size(X)[2]
        pairsMatrix[j,i] = pairsMatrix[i,j] = sort([id2uniprot(ids[j]), id2uniprot(ids[i])])
    end

    mask = ChromNetPaper.upper(M)
    truth = ChromNetPaper.upper(T)[mask]
    scores = ChromNetPaper.upper(X)[mask]
    pairs = ChromNetPaper.upper(pairsMatrix)[mask]

    if numEdges != nothing
        inds = sortperm(-scores)[1:numEdges]
        #println("$numEdges edges chosen.")
        trueInds = inds[truth[inds]]
        return unique(pairs[trueInds])
    elseif scoreThreshold != nothing
        #println(sum(scores .> scoreThreshold), " edges above threshold.")
        trueInds = find(truth .* (scores .> scoreThreshold))
        return unique(pairs[trueInds])
    else
        error("Either numEdges or scoreThreshold must be set!")
    end
end

# function area_under_pr(X::AbstractMatrix, T::Array{Bool,2}, M::Array{Bool,2}; resolution=4000)
#     @assert size(X) == size(T)
#     @assert size(X) == size(M)
#
#     mask = upper(M)
#     scores = upper(X)[mask]
#     truth = upper(T)[mask]
#
#     area_under_pr(truth, scores; resolution=resolution)
# end
#
# function area_under_pr(truth::AbstractVector, predictor::AbstractVector; resolution=4000)
#     rocData = MLBase.roc(round(Int64, truth), float(invperm(sortperm(predictor))), resolution)
#     vals = collect(map(x->(recall(x), -precision(x)), rocData))
#     sort!(vals)
#     xvals = map(x->x[1], vals)
#     yvals = map(x->-x[2], vals)
#     # println(xvals)
#     # println(yvals)
#     area_under_curve([0.0; xvals], [yvals[1]; yvals]) # make sure we extend all the way to zero
# end
# function area_under_curve(x, y) # must be sorted by increasing x
#     area = 0.0
#     lastVal = NaN
#     for i in 2:length(x)
#         v = (y[i-1]+y[i])/2 * (x[i]-x[i-1])
#         if !isnan(v)
#             area += v
#             lastVal = v
#         elseif !isnan(lastVal)
#             area += lastVal
#         end
#     end
#     area
# end
# function area_under_pr(truth::AbstractMatrix, predictor::AbstractMatrix)
#     area_under_pr(abs(upper(truth)) .> 0.01, abs(upper(predictor)))
# end
