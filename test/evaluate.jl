
# upper
X = [1 0.1 0.4 0.2; 0.1 1 0.3 0.5; 0.4 0.3 1 0.6; 0.2 0.5 0.6 1]
@test ChromNetPaper.upper(X) == [0.1, 0.4, 0.3, 0.2, 0.5, 0.6]

# network_enrichment
T = round(Bool, [1 1 1 0; 1 1 0 1; 1 0 1 0; 0 1 0 1])
@test ChromNetPaper.network_enrichment(X, T) == (2/3)/0.5
M = round(Bool, [1 1 1 0; 1 1 1 0; 1 1 1 0; 0 0 0 0])
@test ChromNetPaper.network_enrichment(X, T, M) == (1/2)/(2/3)

# id2uniprot
@test ChromNetPaper.id2uniprot("ENCSR177HDZ") == "P01100"
@test ChromNetPaper.id2uniprot("ENCSR664POU") == "Q04206"
@test ChromNetPaper.id2uniprot("ENCSR459FTB") == "P17480"

# id2celltype
@test ChromNetPaper.id2celltype("ENCSR177HDZ") == "HepG2"

# id2truth
@test ChromNetPaper.id2truth("ENCSR177HDZ", "ENCSR664POU")
@test !ChromNetPaper.id2truth("ENCSR177HDZ", "ENCSR459FTB")

# truth_matrix
ans = [
    true  true  false;
    true  true  false;
    false false true;
]
@test all(ChromNetPaper.truth_matrix(["ENCSR177HDZ", "ENCSR664POU", "ENCSR459FTB"]) .== ans)

# ishistone
@test ChromNetPaper.ishistone("ENCSR449AYM")
@test !ChromNetPaper.ishistone("ENCSR177HDZ")

# area_under_pr
#ChromNetPaper.area_under_pr(X, T, M) # just make sure it runs for now
