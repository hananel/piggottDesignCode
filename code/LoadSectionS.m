function [Cl Cd] = LoadSectionS(alpha, Re, blade)
% general aerodynamic data files are combined arrays of Cl for all reynolds
% and alpha's, and another file for Cd for all reynolds and alpha's
%
% uses interpolation between sections at different radius (in difference
% from LoadSection.m which uses the section in SectionName
% the different section names are at blade, blade.Section at blade.rOrig


% Re color vec
colors = jet(length(Re));

% precautions
alpha = real(alpha);
Re = killNans(Re);
% preallocation
Cl = ones(length(alpha),1);
Cd = ones(length(alpha),1);

for i=1:length(alpha)
    % finding upper bound Section, and lower bound section
    NBot = find(blade.rOrig<=blade.r(i),1,'last');
    rBot = blade.rOrig(NBot);
    NameBot = blade.Section{NBot};
    NTop = find(blade.rOrig>=blade.r(i),1,'first');
    rTop = blade.rOrig(NTop);
    NameTop = blade.Section{NTop};
    SectionPathBot = GetSectionPath(NameBot);
    SectionPathTop = GetSectionPath(NameTop);
    
    % getting section data
    % Bottom
    D = csvread([SectionPathBot 'Cl.txt']); ReDataBot = D(1,2:end);     
    ClDataBot = D(2:end,2:end); alphaDataBot = D(2:end,1); % [deg]
    D = csvread([SectionPathBot 'Cd.txt']); CdDataBot = D(2:end,2:end);
    % Top
    D = csvread([SectionPathTop 'Cl.txt']); ReDataTop = D(1,2:end);     
    ClDataTop = D(2:end,2:end); alphaDataTop = D(2:end,1); % [deg]
    D = csvread([SectionPathTop 'Cd.txt']); CdDataTop = D(2:end,2:end);
    % carrying it all the way - can't make an interpolated matrix yet - not
    % the same range of Re and alpha neccesaeraly...
    
    % precautions
    if Re(i)<max([min(ReDataBot), min(ReDataTop)]) ,Re(i) = max([min(ReDataBot), min(ReDataTop)]); end
    if Re(i)>min([max(ReDataBot), max(ReDataTop)]) ,Re(i) = min([max(ReDataBot), max(ReDataTop)]); end
    
    % interp Cl and Cd
    if alpha(i)>min(max(alphaDataTop),max(alphaDataBot))
        % use NACA0012 data - off point design
        load ../Input/Section/NACA0012/NACA0012.mat
        % shift NACA0012 data to the relevant Reynolds range (LONG)
        % TODO - use linear interpolation of Bot and Top profile, and more reasonable Reynolds choice.
        % at the moment - just taking the Bottom section, with the highest Reynolds. 
        Section.alpha = alphaDataBot; Section.cl = ClDataBot(:,end); Section.cd = CdDataBot(:,end);
        [NACA0012_shift] = shiftNACA0012_byFeature(Section, NACA0012);
        CL = killnans(NACA0012_shift.cl); CD = killnans(NACA0012_shift.cd); Alpha = killnans(NACA0012_shift.alpha);
        
        % other option - make a better fit and save the txt file in the
        % name below
        % D = csvread(['../Input/Section/' SectionName '/' SectionName '2NACA0012.txt']);
        % Alpha = D(:,1); CL = D(:,2); CD = D(:,3);
        
        Cl(i) = interp1(Alpha,CL,alpha(i),'linear');
        Cd(i) = interp1(Alpha,CD,alpha(i),'linear');
    else
        if alpha<max(min(alphaDataBot),min(alphaDataTop)) % no values - extrap (for very high TSR, off design in reality)
            % use highest Reynolds data - simply because i don't have any
            % option to extrap in interp2d ... and this is way off design
            Cl(i) = interp1(alphaDataBot,ClDataBot(:,end),alpha(i),'linear','extrap');
            Cd(i) = interp1(alphaDataBot,CdDataBot(:,end),alpha(i),'linear','extrap');
        else
            % use real data - real point design (TODO add indicator, for validity of data, for future estimate of result)
            if NameTop==NameBot
                Cl(i) = interp2(ReDataBot,alphaDataBot,ClDataBot,Re(i),alpha(i),'linear',0);
                Cd(i) = interp2(ReDataBot,alphaDataBot,CdDataBot,Re(i),alpha(i),'linear',0);
            else
                ClBot = interp2(ReDataBot,alphaDataBot,ClDataBot,Re(i),alpha(i),'linear',0);
                CdBot = interp2(ReDataBot,alphaDataBot,CdDataBot,Re(i),alpha(i),'linear',0);
                ClTop = interp2(ReDataTop,alphaDataTop,ClDataTop,Re(i),alpha(i),'linear',0);
                CdTop = interp2(ReDataTop,alphaDataTop,CdDataTop,Re(i),alpha(i),'linear',0);
                Cl(i) = ClBot + (blade.r(i)-rBot)/(rTop-rBot)*(ClTop - ClBot);
                Cd(i) = CdBot + (blade.r(i)-rBot)/(rTop-rBot)*(CdTop - CdBot);
            end
        end
    end
end

%%
% Reynolds interpolation for each section - producing number
% interpolating between sections - on this number and rOrig vs. r
% if the angle is too high, doing as above, only for all alphas - to produce Section.alph etc. for NACA0012_shift
