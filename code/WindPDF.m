function windPdfInRange = WindPDF(height, windRange,a,H0,WBLscale,WBLshape)
% yearlyProduction = 365*24*dw*sum(elecPower.*windPdfInRange)
% the input here is the windspeed  values in reference to hub height
% height is the tower height
% H0 = measurement height
% WBLSCALE = Weibull shape parameter
% ALPHA = Blausius coefficient; .17

K = (H0/height)^a; 			
windRangeH0 = windRange*K; 		
windPdfInRange = wblpdf(windRangeH0, WBLscale, WBLshape)*K;
Uavg = WBLscale*gamma(1+1/WBLscale)