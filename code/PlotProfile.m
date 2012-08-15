function [x,y,name] = PlotProfile(fileName)
% fileName is a txt file containing a header which is the profile name, and 2
% columns, of x and y coordinates respectfully
% output: 
% name: section name
% x: section x data
% y: section y data

M = dlmread(fileName,' ',1,0);
% killing zeroes from mis-reading by dlmread
for i=1:size(M,1)
    for j=2:size(M,2)
        if or(M(i,j)>0,M(i,j)<0)
            M(i,2) = M(i,j);
            break;
        end
    end
end
x = M(:,1); y = M(:,2);
fid = fopen(fileName);
name = fgets(fid);

plot(x,y); axis equal
title(name);

