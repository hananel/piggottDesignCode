function minimum = test_NewBlade(inp)
% running OptimizeNewBlade
tic
inp_CometME42Windy
[minimum] = NewBlade(inp);
minutes = toc/60