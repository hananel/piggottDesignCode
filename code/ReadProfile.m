function [x,y,name] = ReadProfile(fileName)
% fileName is a txt file containing a header which is the profile name, and 2
% columns, of x and y coordinates respectfully
% output:
% name: section name
% x: section x data
% y: section y data

M = dlmread(fileName,' ',1,0);
if M(:,2)==0
    M = dlmread(fileName,'',0,0);
    name = fileName;
else
    fid = fopen(fileName);
    name = fgets(fid);
end
x = M(:,1); y = M(:,2);
plot(x,y); axis equal
title(name);

