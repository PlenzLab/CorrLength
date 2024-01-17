function CorrelationLengthAnalysis(timeseries, coords, varargin)
% CorrelationLengthAnalysis(timeseries, coords, [optional args])
%
% This function performs spatial correlation analysis as shown in
% Ribeiro et al. 2024 Cell Reports paper titled "Trial-by-trial
% variability in cortical responses exhibits scaling in spatial 
% correlations predicted from critical dynamics"
%
% Written by Tiago Lins Ribeiro (tiagolinsribeiro@gmail.com) 2024-01-17
%
% INPUT:
% timeseries  Npoints x Nunits matrix with the time series of each unit
% coords      Nunits x 2 matrix with the location (x, y) of each unit
%
% OPTIONAL INPUT:
% -d  dsfactor     downsampling factor multiplier (default 1)
% -nL NLs          number of sampling sizes employed (default 15)
% -wf wf           fraction of the max # of windows used (default 0.05)
% -nb Nbins        number of bins in the distance domain (default 20)
% -L  maxL         max window size (default 0, calculates from the data)
% -l  minL         min window size (default 0, 10% of maxL)
% -o  outfile      name (and path) of the output file (default CorrL)
%
% Example usage:
% CorrelationLengthAnalysis(data, coords, '-NLs', 20, '-wf', 0.01, '-d', 2)


% =========================================================================
% Set up variables for the code to run.
% =========================================================================

% Exit function if the required arguments are not provided
if nargin < 2, error('Insuficient input files provided.'); end

% Default values (see explanation above)
dsfactor = 1;
NLs      = 15;
wf       = 0.05;
Nbins    = 20;
maxL     = 0;
minL     = 0;
outfile  = 'CorrL';

% Parse varargin to read all optional inputs provided
while ~isempty(varargin)
    switch varargin{1}
         
        case '-d'
            dsfactor = varargin{2};
            varargin(1:2) = [];
        
        case '-nL'
            NLs = varargin{2};
            varargin(1:2) = [];
        
        case '-wf'
            wf = varargin{2};
            varargin(1:2) = [];
        
        case '-nb'
            Nbins = varargin{2};
            varargin(1:2) = [];
        
        case '-L'
            maxL = varargin{2};
            varargin(1:2) = [];
        
        case '-l'
            minL = varargin{2};
            varargin(1:2) = [];
        
        case '-o'
            outfile = varargin{2};
            varargin(1:2) = [];
        
        otherwise
            error(['Unexpected option: ' varargin{1}])
    end
end

% Create the output file if it doesn't exist
if ~exist(outfile,'file')
    createdon = strrep(datestr(now),' ','-'); 
    save(outfile,'createdon', 'NLs', 'Nbins', 'wf', 'dsfactor','-v7.3');
end

% =========================================================================
% Process the inputs and organizes into new variables. 
% =========================================================================

xc      = coords(:,1);
yc      = coords(:,2);

% Make sure window starts at zero
xc = xc - min(xc); 
yc = yc - min(yc);

Npoints = size(timeseries, 1);
Nunits  = size(timeseries, 2);
   
% Downsample the time series
Npoints_ds = ( Npoints + mod(-Npoints,dsfactor) )/dsfactor;
ds_series  = zeros(Npoints_ds, Nunits);
for u = 1:Nunits
    ds_series(:,u) = nanmean(reshape([timeseries(:,u);...
        nan(mod(-Npoints,dsfactor),1)],dsfactor,[]));
end

save(outfile,...
	'xc',...
	'yc',...
	'ds_series',...
	'Npoints',...
	'Nunits',...
	'-append');

% =========================================================================
% Define parameters of the spatial analysis.
% =========================================================================

% Obtain sample window sizes with logarithmic spacing from minL to maxL
L     = maxL; if L == 0, L = max(coords(:)); end
Lsmin = minL; if Lsmin == 0, Lsmin = 0.1*L; end
Lsmax = L;
dL    = (Lsmax/Lsmin)^(1/(NLs-1));
Ls    = 0:NLs-1;
Ls    = round(Lsmin*dL.^Ls);
Ls    = unique(Ls);
NLs   = length(Ls);

% Obtain specific windows to be sampled by randomly picking from the pool
windex   = cell(NLs, 1);
Nwindows = zeros(NLs, 1);
for s = 1:NLs
    maxNw       = (L + 1 - Ls(s))^2; % see Note1 below
    Nwindows(s) = ceil(wf*maxNw); % see Note2 below
    windex{s}   = randperm(maxNw, Nwindows(s));
end
% Note1: step of 1 when calculating maximum number of possible windows
% might be incompatible with data. Adjust accordingly.
% Note2: make sure the selected fraction of sampled windows results in above
% reasonable number of samples in terms of processing time and statistics.
% Might have to adjust based on spatial distribution of units. 

% Save variables to the output file
save(outfile,...
    'L',...
    'NLs',...
    'Ls',...
    'Nwindows',...
    'windex',...
    '-append');

% =========================================================================
% Obtain distance and correlation for unit pairs for each sampled window.
% =========================================================================

D.s = cell(NLs, 1);
C.s = cell(NLs, 1);
for s = 1:NLs
    D.s{s}.b = cell(Nbins, 1);
    C.s{s}.b = cell(Nbins, 1);
    for b = 1:Nbins
        D.s{s}.b{b}.w = cell(Nwindows(s), 1);
        C.s{s}.b{b}.w = cell(Nwindows(s), 1);
    end
	
	for w = 1:Nwindows(s)
	
		% Excludes units outside the sampling window
		xmin     = mod(windex{s}(w)-1, L + 1 - Ls(s)) + 1;
		xmax     = xmin + Ls(s) - 1;
		ymin     = floor((windex{s}(w)-1)/(L+1-Ls(s))) + 1;
		ymax     = ymin + Ls(s) - 1;
		sample   = find(xc' >= xmin & xc' <= xmax & yc' >= ymin & yc' <= ymax);
		Nsampled = length(sample);
	
		% !Proceed only if there are at least 4 units in the sampling window!
		if Nsampled >= 4
			
			% Distance matrix
			xs = xc(sample);
			ys = yc(sample);
			Dm = zeros(Nsampled, Nsampled);
			for u1 = 1:Nsampled
				for u2 = 1:Nsampled
					Dm(u1,u2) = sqrt((xs(u1)-xs(u2))^2 + (ys(u1)-ys(u2))^2);
				end
			end			
			
			% Obtain the correlation matrix after population mean subtraction
			avgS = mean(ds_series(1:end-1,sample), 2);
			S    = ds_series(1:end-1,sample) - avgS;
			R    = corrcoef(S);
			
			% Compute covariance versus distance
			d  = cell(Nbins, 1);
			c  = cell(Nbins, 1);
			R  = R(:);
			Dm = Dm(:);
			dD = (Ls(s) - 0.5)/Nbins;
			for b = 1:Nbins
				d{b} = [d{b} Dm(Dm > (b-1)*dD+0.5 & Dm <= b*dD+0.5)];
				c{b} = [c{b} R(Dm > (b-1)*dD+0.5 & Dm <= b*dD+0.5)];
			end
			
			for b = 1:Nbins
				D.s{s}.b{b}.w{w} = d{b};
				C.s{s}.b{b}.w{w} = c{b};
			end
			clear Dm S R d c;

		end
	
	end
end

% =========================================================================
% Obtain average correlation functions and correlation length.
% =========================================================================

x0   = linspace(0, L, Nbins);
y0   = zeros(Nbins, 1);
xdat = cell(NLs, 1);
ydat = cell(NLs, 1);
fo   = fitoptions(...
    'Method',     'NonlinearLeastSquares',...
    'Lower',      [-1 0 0 -100],...
    'Upper',      [1 100 100 100],...
    'StartPoint', [-0.1 1 1 0]);
ft   = fittype('a + b*exp(-c*x) - d*x', 'options', fo); % see Note3 below
% Note3: this function and respective parameters are used to fit the
% correlation function and get a better estimate of its zero-crossing.
% Particular shape and parameters will depend on data, please check.

avgD = cell(NLs, 1);
avgC = cell(NLs, 1);
stdC = cell(NLs, 1);
semC = cell(NLs, 1);
CL   = nan(NLs, 1);
cfit = cell(NLs, 1);
r0   = nan(NLs, 1);
for s = 1:NLs
    
    xdat{s} = linspace(0, Ls(s), 10000);
    avgD{s} = cell(Nbins, 1);
    avgC{s} = cell(Nbins, 1);
    stdC{s} = cell(Nbins, 1);
    semC{s} = cell(Nbins, 1);
    
    for b = 1:Nbins
        
        for w = 1:Nwindows(s)
        
            % Pool data from all the windows together
            avgD{s}{b} = [avgD{s}{b}; D.s{s}.b{b}.w{w}];
            avgC{s}{b} = [avgC{s}{b}; C.s{s}.b{b}.w{w}];
        
        end
        
        % Average data from all windows
        avgD{s}{b} = nanmean(avgD{s}{b});
		stdC{s}{b} = nanstd(avgC{s}{b},1);
		semC{s}{b} = stdC{s}{b}/sqrt(length(avgC{s}{b}));
        avgC{s}{b} = nanmean(avgC{s}{b});

        
    end
    
    % Reformat data
    avgD{s} = cell2mat(avgD{s});
    stdC{s} = cell2mat(stdC{s});
    semC{s} = cell2mat(semC{s});
    avgC{s} = cell2mat(avgC{s});
    
	
	% Calculate the zero-crossing of correlation for this window size
    xfit = avgD{s}(~isnan(avgC{s}));
    yfit = avgC{s}(~isnan(avgC{s}));
    if length(xfit(xfit<=Ls(s))) >= 4
        cfit{s} = fit(xfit,yfit,ft,'Exclude',xfit>Ls(s));
        ydat{s} = cfit{s}(xdat{s});
        [cross0, ~] = intersections(xdat{s}, ydat{s}, x0, y0);
    else
        cross0 = [];
    end
    if ~isempty(cross0), r0(s) = cross0(1); end

    % Calculate the correlation length for this window size
	CL(s) = sqrt(nansum(avgC{s}(avgC{s}>0).*(avgD{s}(avgC{s}>0).^2))/nansum(avgC{s}(avgC{s}>0))); % see Note4 below
	% Note4: Depending on how noisy the correlation functions are, it might
	% be better to obtain the correlation length from the correlation function
	% fits, such as:
	% CL(s) = sqrt(nansum(ydat{s}(xdat{s}<=r0(s)).*(xdat{s}(xdat{s}<=r0(s)).^2))/nansum(xdat{s}(xdat{s}<=r0(s))));


end


% Append the variables to the file provided
save(outfile,...
    'avgD',...
    'avgC',...
    'stdC',...
    'semC',...
    'xdat',...
    'ydat',...
    'cfit',...
    'CL',...
    'r0',...
    '-append');


% =========================================================================
% Perform collapse analysis.
% =========================================================================

% Calculate the slope at zero crossing and the exponent gamma
dCdD0 = nan(NLs, 1);
for s = 1:NLs
    if ~isempty(cfit{s})
        dCdD0(s) = abs(differentiate(cfit{s},r0(s)));
    end
end
if length(r0(~isnan(r0))) >= 4
    plfit = fit(r0(~isnan(r0)),dCdD0(~isnan(r0)).*r0(~isnan(r0)),'a*x^b');
    gamma = -1*plfit.b;
else
    gamma = 0;
end

% Obtain the rescaled x axis
sclD = zeros(NLs, Nbins);
for s = 1:NLs
    for b = 1:Nbins
        sclD(s,b) = avgD{s}(b)/r0(s);
    end
end

% Obtain the rescaled y axis
sclC = zeros(NLs,Nbins);
for s = 1:NLs
    for b = 1:Nbins
        sclC(s,b) = avgC{s}(b)*(r0(s).^gamma);
    end
end

% Obtain average and s.d. rescaled correlations
dbin = 2/Nbins;
sclDbin = zeros(Nbins, NLs);
sclCbin = zeros(Nbins, NLs);
for s = 1:NLs
    for b = 1:Nbins
        sclDbin(b,s) = nanmean(sclD(s, sclD(s,:)>=(b-1)*dbin & ...
            sclD(s,:)<b*dbin));
        sclCbin(b,s) = nanmean(sclC(s, sclD(s,:)>=(b-1)*dbin & ...
            sclD(s,:)<b*dbin));
    end
end
sclDavg = nanmean(sclDbin,2);
sclCavg = nanmean(sclCbin,2);
sclCstd = nanstd(sclCbin,1,2);


% Calculate collapse quality
dist   = zeros(NLs, Nbins);
for s = 1:NLs
    for b = 1:Nbins
        dist(s,b)   = (sclCbin(b,s) - sclCavg(b))/sclCavg(b);
    end
end
collapseQ = 100*nanmean(abs(nanmean(dist(sclD<1),1)));

% Calculate area under the curve
collapseA = zeros(NLs,1);
for s = 1:NLs
    if ~isempty(cfit{s})
        collapseA(s) = integrate(cfit{s},r0(s),0); % Susceptibility
    end
end


% Append the variables to the file provided
save(outfile,...
    'sclD',...
    'sclC',...
    'sclDbin',...
    'sclCbin',...
    'sclDavg',...
    'sclCavg',...
    'sclCstd',...
    'dCdD0',...
    'plfit',...
    'gamma',...
    'collapseQ',...
    'collapseA',...
    '-append');

% =========================================================================
% Prepare figures.
% =========================================================================

figure;
cmap = colormap(copper(NLs));

% Plot the average + s.e.m. connected correlation as function of distance
% for different sampling window sizes
subplot(1,3,1);
xmin = 0;
xmax = 0.8*L;
for s = 1:NLs
    errorbar(avgD{s}(~isnan(avgD{s})),...
        avgC{s}(~isnan(avgD{s})), ...
        semC{s}(~isnan(avgD{s})),...
        '-o','Linewidth',2.5,'MarkerSize',6.,...
        'Color',cmap(s,:));
    hold on;
    ydat = cfit{s}(xdat{s});
    plot(xdat{s}, ydat,...
        '--','Linewidth',2.5,'MarkerSize',6.,...
        'Color',cmap(s,:));
end
xlim([xmin xmax]);
ylim([-0.25 1]);
xlabel('Distance');
ylabel('Correlation');
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 2;
pbaspect([1 1 1]);

% Plot the rescaled correlation as function of distance
% for different sampling window sizes
subplot(1,3,2);
xmin = 0;
xmax = 2;
for s = 1:NLs
    plot(sclD(s,~isnan(sclD(s,:))),...
        sclC(s,~isnan(sclD(s,:))), ...
        'o','MarkerSize',6.,'Color',cmap(s,:));
    hold on;
end
errorbar(sclDavg, sclCavg, sclCstd,...
    '-k','Linewidth',3.5);
xlim([xmin xmax]);
xlabel('Distance/r0');
ylabel('Correlation*r0^g');
hline = refline([0 0]);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 2;
pbaspect([1 1 1]);


% Plot the correlation length as function of sampling window size
subplot(1,3,3);
plot(Ls, CL, 'o', 'MarkerSize', 6., 'Color', cmap(NLs,:)); hold on;
plot(Ls, 0.45*Ls, '--k');
xlabel('Window size');
ylabel('\xi');
pbaspect([1 1 1]);

saveas(gcf, [outfile '.fig']);
saveas(gcf, [outfile '.tif']);

end

