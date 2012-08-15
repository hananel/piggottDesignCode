function [CpData CtData TSR VData] = loadCpCtData(inp)

CpData = csvread(['../' inp.CpCt_fileDirectory inp.Cp_filename]);
CtData = csvread(['../' inp.CpCt_fileDirectory inp.Ct_filename]);
VData = CpData(1,2:end);
TSR = CpData(2:end,1);
LegendText = arrayfun(@(n) {sprintf('%d', n)}, (VData));

colorData = 'rgbmckrgbmck';
if nargin==0
    figure(4); hold on;
    for i=1:length(VData)
        subplot(211); hold on;
        plot(TSR,CpData(2:end,1+i),colorData(i));
        subplot(212); hold on;
        plot(TSR,CtData(2:end,1+i),colorData(i));
    end
    subplot(211); ylabel('Cp'); xlabel('TSR'); title('WindyLight03A18 optimal TSR2 2/1/11')
    axis([0 max(TSR) 0 0.593])
    subplot(212); ylabel('Ct'); xlabel('TSR');
    legend(LegendText,'Location','NorthWest')
end