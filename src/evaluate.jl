
export network_enrichment, id2uniprot, id2truth, id2celltype, ishistone, truth_matrix

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

id2uniprotDict = Dict()
"Get the uniprot id of a given experiment id."
function id2uniprot(id)
    global id2uniprotDict

    # lazy load the data
    if length(id2uniprotDict) == 0
        f = open(joinpath(globalDataDir, "id2uniprot.txt"))
        for line in eachline(f)
            parts = split(strip(line), '\t')
            id2uniprotDict[parts[1]] = parts[2]
        end
        close(f)
    end

    get(id2uniprotDict, id, "")
end


id2cellTypeDict = Dict()
"Get the cell type of a given experiment id."
function id2celltype(id)
    global id2cellTypeDict

    # lazy load the data
    if length(id2cellTypeDict) == 0
        f = open(joinpath(globalDataDir, "id2celltype.txt"))
        for line in eachline(f)
            parts = split(strip(line), '\t')
            id2cellTypeDict[parts[1]] = parts[2]
        end
        close(f)
    end

    get(id2cellTypeDict, id, "")
end
