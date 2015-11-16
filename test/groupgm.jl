using ChromNet

# project_groupgm
IC = [1 0.1 0.4 0.2; 0.1 1 0.3 0.5; 0.4 0.3 1 0.6; 0.2 0.5 0.6 1]
header = ["A1", "A2", "B1", "B2"]
groups = build_groups(inv(IC), header)
G,headerG = build_groupgm(IC, header, groups)

ans = [
    1.0 0.9 0.6 0.6;
    0.9 1.0 0.9 0.9;
    0.6 0.9 1.0 0.6;
    0.6 0.9 0.6 1.0
]
@test all((ChromNetPaper.project_groupgm(G, header, groups, x->true) .- ans) .< 1e-8)

ans = [
    1.0 0.1 0.6 0.6;
    0.1 1.0 0.8 0.8;
    0.6 0.8 1.0 0.6;
    0.6 0.8 0.6 1.0
]
@test all((ChromNetPaper.project_groupgm(G, header, groups, x->x[1] < 0.1) .- ans) .< 1e-8)

# run with default group filter
ChromNetPaper.project_groupgm(G, header, groups)
