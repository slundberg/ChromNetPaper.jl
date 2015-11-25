using ChromNet

export project_groupgm, edge_groupgm

"""
This projects a group graphical model down to a network among individual variables.

Projection is necessary in order to compare with standard modeling methods that
only connect individual datasets. The group filter determines which groups should
be considered for the projection (this may exclude large groups), and then the
edge merge function determines how to merge two edge values.
"""
function project_groupgm(G::AbstractMatrix, header::AbstractArray, groups::AbstractArray, groupFilter=nothing, edgeMerge::Function=(x,y)->abs(x) > abs(y) ? x : y)

    # by default only allow groups targeting the same factor
    function same_factor_filter(x)
        parts = split(x[2], ' ')
        referenceUid = ChromNetPaper.id2uniprot(parts[1])
        for i in 2:length(parts)
            if referenceUid != ChromNetPaper.id2uniprot(parts[i])
                return false
            end
        end
        true
    end
    if groupFilter == nothing
        groupFilter = same_factor_filter
    end

    # create a map from header ids to indexes
    indexMap = Dict{ASCIIString,Int64}()
    for i in 1:length(header)
        indexMap[header[i]] = i
    end

    # fill in our projected matrix
    X = NaN*zeros(length(header),length(header))
    indexGroups = [[indexMap[x] for x in split(g[2])] for g in groups]
    for i in 1:length(indexGroups)
        if !groupFilter(groups[i]) continue end

        # fill in diagonal
        if length(indexGroups[i]) == 1
            ind = indexGroups[i][1]
            X[ind,ind] = G[i,i]
        end

        for j in i+1:length(indexGroups)
            if !groupFilter(groups[j]) || issubset(indexGroups[i], indexGroups[j]) || issubset(indexGroups[j], indexGroups[i])
                continue
            end

            for gi in indexGroups[i], gj in indexGroups[j]
                X[gj,gi] = X[gi,gj] = isnan(X[gj,gi]) ? G[j,i] : edgeMerge(X[gj,gi], G[j,i])
            end
        end
    end

    X
end

function edge_groupgm(C, ids, invFunc=inv; checkGroup=nothing)
    groups = build_groups(C, ids)
    G,headerG = build_groupgm(invFunc(C), ids, groups)
    project_groupgm(G, ids, groups, checkGroup)
end
