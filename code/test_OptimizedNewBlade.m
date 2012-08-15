function minimum = test_OptimizedNewBlade(inp)
% running OptimizeNewBlade
tic
inp_Windylight_normal_fucking_turbine
[minimum] = OptimizeNewBlade(2.5,inp);
minutes = toc/60