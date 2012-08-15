function plotXfoil(color)

% pick xFoil txt file
[filename, pathname, filterindex] = uigetfile({'*.txt';'*.*'}, ...
    'Pick a xFoil output file');

d = dlmread([pathname filename],'',12,0);
a = d(:,1); cl = d(:,2); cd = d(:,3);

subplot(121); hold on; plot(a,cl,color,a,cd,color)
xlabel('alpha'); ylabel('cl, cd');
subplot(122); hold on; plot(a,cl./cd,color)
xlabel('alpha'); ylabel('cl/cd')