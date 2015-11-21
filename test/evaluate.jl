
# upper
X = [1 0.1 0.4 0.2; 0.1 1 0.3 0.5; 0.4 0.3 1 0.6; 0.2 0.5 0.6 1]
@test ChromNetPaper.upper(X) == [0.1, 0.4, 0.3, 0.2, 0.5, 0.6]

# network_enrichment
T = round(Bool, [1 1 1 0; 1 1 0 1; 1 0 1 0; 0 1 0 1])
@test network_enrichment(X, T) == (2/3)/0.5
M = round(Bool, [1 1 1 0; 1 1 1 0; 1 1 1 0; 0 0 0 0])
@test network_enrichment(X, T, M) == (1/2)/(2/3)

# id2uniprot
@test id2uniprot("ENCSR177HDZ") == "P01100"
@test id2uniprot("ENCSR664POU") == "Q04206"
@test id2uniprot("ENCSR459FTB") == "P17480"

# id2celltype
@test id2celltype("ENCSR177HDZ") == "HepG2"

# id2target
@test id2target("ENCSR177HDZ") == "FOS"

# id2treatments
@test id2treatments("ENCSR177HDZ") == "None"

# id2truth
@test id2truth("ENCSR177HDZ", "ENCSR664POU")
@test !id2truth("ENCSR177HDZ", "ENCSR459FTB")

# truth_matrix
ans = [
    true  true  false;
    true  true  false;
    false false true;
]
@test all(truth_matrix(["ENCSR177HDZ", "ENCSR664POU", "ENCSR459FTB"]) .== ans)

# ishistone
@test ishistone("ENCSR449AYM")
@test !ishistone("ENCSR177HDZ")

# mask_matrix
ans = [
    false false false;
    false false true;
    false true  false;
]
@test all(mask_matrix("within_all", ["ENCSR177HDZ", "ENCSR664POU", "ENCSR459FTB"]) .== ans)

# enrichment_rank
mask = ChromNetPaper.upper(M)
scores = ChromNetPaper.upper(X)[mask]
truth = ChromNetPaper.upper(T)[mask]
x,y = enrichment_rank(truth, scores)
@test maximum(abs(x .- [1,2,3])) < 1e-8
@test maximum(abs(y .- [(1/1)/(2/3), (1/2)/(2/3), (2/3)/(2/3)])) < 1e-8

# network_enrichment_density
x2,y2 = network_enrichment_density(X, T, M)
@test maximum(abs(y2 .- y)) < 1e-8
@test abs(x2[end] - 1.0) < 1e-8

# network_enrichment_rank
x2,y2 = network_enrichment_rank(X, T, M)
@test maximum(abs(y2 .- y)) < 1e-8
@test abs(x2[end] - length(x2)) < 1e-8
