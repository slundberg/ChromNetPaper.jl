
export project_groupgm

"""
This projects a group graphical model down to a network among individual variables.

Projection is nessecary in order to compare with standard modeling methods that
only connect individual datasets. The group filter determines which groups should
be considered for the projection (this may exclude large groups), and then the
edge merge function determines how to merge two edge values.
"""
function project_groupgm(G::AbstractMatrix, header::AbstractArray, groups::AbstractArray, groupFilter::Function, edgeMerge::Function=(x,y)->abs(x) > abs(y) ? x : y)

    # create a map from header ids to indexes
    indexMap = Dict{ASCIIString,Int64}()
    for i in 1:length(header)
        indexMap[header[i]] = i
    end

    # fill in our projected matrix
    X = NaN*zeros(length(header),length(header))
    indexGroups = collect(map(g->collect(map(x->indexMap[x], split(g[2]))), groups))
    for i in 1:length(indexGroups)
        if !groupFilter(groups[i]) continue end

        # fill in diagonal
        if length(indexGroups[i]) == 1
            ind = indexGroups[i]
            X[ind,ind] = G[i,i]
        end

        for j in i:length(indexGroups)
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
