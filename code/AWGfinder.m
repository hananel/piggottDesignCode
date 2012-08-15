function [AWG,sqmm,AWG_sqmm,wireN] = AWGfinder(area);
% takes a sqmm size (size) and findes the nearest smaller AWG size from table

Table = [ 40  0.0787  0.5189
          39  0.0889  0.4116
          38  0.1016  0.3269
          37  0.1143  0.2588
          36  0.1270  0.2047
          35  0.1422  0.1624
          34  0.1600  0.1281
          33  0.1803  0.1022
          32  0.2032  0.0804
          31  0.2261  0.0647
          30  0.2540  0.0507
          29  0.2870  0.0401
          28  0.3200  0.0324
          27  0.3607  0.0255
          26  0.4039  0.0201
          25  0.4547  0.0159
          24  0.5105  0.0127
          23  0.5740  0.0103
          22  0.6452  0.0081
          21  0.7239  0.0062
          20  0.8128  0.0049                    
          19  0.91    0.65
          18  1.02    0.82
          17  1.15    1.04
          16  1.29    1.31
          15  1.45    1.65
          14  1.63    2.08
          13  1.83    2.62
          12  2.05    3.31
          11  2.31    4.17
          10  2.59    5.26
          9   2.91    6.62
          8   3.26    6.62
          7   3.67    10.6
          6   4.11    13.3
          5   4.62    16.8
          4   5.19    21.2
          3   5.83    26.7
          2   6.54    33.6
          1   7.35    42.4];
% correcting for enamel
Table(:,2) = Table(:,2)+0.06;
Table(:,3) = pi*Table(:,2).^2/4;

% searching for the fit
% assuming - if the result is bigger then AWG#14 (thath is - smaller AWG)
% then using two wires in parralel
wireN = 1;
AWG = Table(find(Table(:,3)<area,1,'last'),1);
sqmm = Table(find(Table(:,1) == AWG,1,'first'),3);
if AWG<10
    wireN = 2;
    AWG = Table(find(Table(:,3) < area/2,1,'last'),1);
    sqmm = Table(find(Table(:,1) == AWG,1,'first'),3)*2;
end
AWG_sqmm = sqmm/area;