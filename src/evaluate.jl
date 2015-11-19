using MLBase
using JSON

export
    network_enrichment,
    id2uniprot,
    id2truth,
    id2celltype,
    id2treatments,
    id2target,
    ishistone,
    truth_matrix,
    mask_matrix,
    unique_ppi_pairs,
    apply_to_celltypes,
    network_enrichment_density,
    enrichment_rank

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
function network_enrichment(X::AbstractMatrix, T::Array{Bool,2}, M::Array{Bool,2}; numEdges=:num_true)
    @assert size(X) == size(T)
    @assert size(X) == size(M)

    mask = upper(M)
    scores = upper(X)[mask]
    truth = upper(T)[mask]
    if numEdges == :num_true
        numEdges = sum(truth)
    end

    ratioCorrect = sum(truth[sortperm(-scores)][1:numEdges])/numEdges
    ratioCorrect / (sum(truth)/length(truth))
end
function network_enrichment(X::AbstractMatrix, T::Array{Bool,2}; numEdges=:num_true)
    network_enrichment(X, T, ones(Bool, size(X)...); numEdges=numEdges)
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
function enrichment_rank(truth, pred)
    enrichment_rank(truth[sortperm(pred, rev=true)])
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

function network_enrichment_density(X::AbstractMatrix, T::Array{Bool,2}, M::Array{Bool,2})
    @assert size(X) == size(T)
    @assert size(X) == size(M)

    mask = upper(M)
    scores = upper(X)[mask]
    truth = upper(T)[mask]
    x,y = enrichment_rank(truth, scores)
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
    global id2truthDict

    if (id2uniprot(id1) == id2uniprot(id2)) && (id2uniprot(id1) != "")
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

    get(id2truthDict, (id2uniprot(id1), id2uniprot(id2)), false)
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

function mask_matrix(maskType, ids; excludeCellType=nothing, includeCrossEdges=false)
    P = length(ids)
    M = zeros(Bool, P, P)
    for i in 1:P, j in i:P

        # exclude connections involving histones
        exclude = ishistone(ids[i]) || ishistone(ids[j])

        if excludeCellType != nothing
            exclude = exclude || id2celltype(ids[i]) == excludeCellType || id2celltype(ids[j]) == excludeCellType
        end

        # exclude same uniprot target connections
        exclude = exclude || (id2uniprot(ids[i]) == id2uniprot(ids[j]))

        # restrict by cellType
        if maskType == "within_all" || maskType == "within_all_separate"
            exclude = exclude || id2celltype(ids[i]) != id2celltype(ids[j])
        elseif maskType == "between_all"
            exclude = exclude || id2celltype(ids[i]) == id2celltype(ids[j])
        elseif maskType == "all"
            # no filter
        else
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
        println("$numEdges edges chosen.")
        trueInds = inds[truth[inds]]
        return unique(pairs[trueInds])
    elseif scoreThreshold != nothing
        println(sum(scores .> scoreThreshold), " edges above threshold.")
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
