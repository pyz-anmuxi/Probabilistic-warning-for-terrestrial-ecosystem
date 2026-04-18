function [D1, D2, D4, ErrorD1, ErrorD2, N, C] = LangevinReconst(data, L, R, bins, Tau, dt, method, title, useparfor)
%res = LangevinReconst(data, L, R, bins, Tau, dt, method, title, useparfor)
%Implements Langevin Reconstruction
%
%Inputs:
%   data      the data set (evenly sampled time series) (can also take
%             segmented data (see makesegments)
%   L and R   are, respectively, the lower and upper limit considered for the states being spanned by the range
%             of data. Normally, L=min(data) and R=max(data) are good choices but for limited data one might consider a slightly
%             bigger value for L and a slightly smaller value for R.
%   bins      number of bins for the reconstruction Should be carefully chosen depending on the size of data. This number should not be so big so that
%             there would be little data in each bin and it should not be so small so that we are left with few bins. In the following citation, it is recommended to choose
%             bins in such a way that each bin will, at least, contain 100 observations).
%   Tau       a list of time lags considered (integer values) (we chose 1:5. For high resolution data we can consider bigger number of time lags).
%   dt        time step of the data series
%   method    is either 'Nadaraya-Watson' or empty, 'Nadaraya-Watson' is recommended
%             as the result is smoothed.  £¨nonparametric estimation£©  
%
%Outputs:
%   D1        vector of estimated drift function for each bin (first Kramers-Moyal coefficient)
%   D2        vector of estimated diffusion function for each bin (second Kramers-Moyal coefficient)
%   D4        vector of estimated fourth Kramers-Moyal coefficient for each bin (should be small)
%   ErrorD1   vector of estimated error in D1 expressed as standard deviation for each bin.
%   ErrorD2   vector of estimated error in D2 expressed as standard deviation for each bin.
%   N         vector representing the number of data per bin.
%   C         vector of bin centers (So, if you want to plot drift or diffusion functions then the the proper commands
%             are plot(C,D1,'.k') and plot(C,D2,'.k')), respectively.
%
%   this function requires the Curve Fit Toolbox to be installed
%
%   Implemented in Matlab by Babak M.S. Ariani
%   Latest version: https://git.wageningenur.nl/sparcs/langevinreconstruct
%   You can alternatively use the R package with instructions about the method in
%     Rinn, P., Lind, P. G., Whter, M. & Peinke, J. The Langevin Approach: An R Package for Modeling Markov Processes. arXiv preprint arXiv:1603.02036 (2016).
if nargin < 1
    error('LangevinRecostr:nodata', 'No data specified');
end
if nargin < 2 || isempty(L)
    L = min(data);
end
if nargin < 3 || isempty(R)
    R = max(data);
end
if nargin < 4
    bins = 100;
end
if nargin < 5
    Tau = 1:5;
end
if nargin < 6
    dt = 1;
end
if nargin < 7
    method = 'Nadaraya-Watson';
end
if nargin < 8
    title = '';
end
if nargin < 9
    useparfor = true;
end
if useparfor
    numWorkers = Inf;  %inf: ÎÞÇî´ó
else
    numWorkers = 1;
end

if size(data, 2) < size(data, 1)
    data = data.';   %1*60000
end
% ** this is an important difference with some other implementations
% Before selecting the data to be used, we make segments with the history
% for each data point. These segments save the near history of each data point
% per row:
% [x(t) x(t-1), x(t-2),...x(t-(max(Tau)+1))]
% (if one of these values do not exist they are replaced by NAN
% After this we make a selection of the data points (where the
% non-selected data can still be used for the tau.
% Outliers should be replaced by NAN on beforehand.
% We can also use these segments for block bootstrapping
if size(data, 1) == 1
    data = makesegments(data, max(Tau));
end

dat = data(:, data(1, :) >= L & data(1, :) <= R); % here we limit the data to lower and upper bounds L and R
% while keeping the autocorrelation structure correct
x = linspace(L, R, bins + 1);
dx = x(2) - x(1);
C = x + dx / 2;
C = C(1:end - 1); % C is the bin centers.
[N, ~] = histcounts(dat(1, :), x);

% Estimating conditional moments
M1 = zeros(1, bins);
M2 = M1;
M4 = M1;
dat(:, dat(1, :) == R) = R - 10^(-6); % By this trick all data having the maximum value of R will be included in the last bin!
A1 = zeros(length(Tau), bins);
A2 = A1;
A4 = A1; % Ai Matrix contain the Mi vectors for different values of Tau.
ErrorM1 = zeros(length(Tau), bins);
ErrorM2 = ErrorM1; % Error1 and Error2 are variance of estimmated M1 and M2 coefficients for different lags and bins.

% This is the direct way to estimate conditional methods.
if strcmp(method, 'Nadaraya-Watson') == 0
    for tau = Tau
        for i = 1:bins % Here, I use Matlab parallel computing toolbox. If you do not have it just use 'for loop' not 'parfor loop'
            dat1 = dat; % broaccast variables can have overhead
            x1 = x;
            I = dat1(1, :) >= x1(i) & dat1(1, :) < x1(i + 1);
            B = dat1(tau + 1, I) - dat1(1, I);
            ndx = ~isnan(B);
            B = B(ndx);
            M1(i) = mean(B); % M1 = <x(t+tau)-x(t)>|x(t)=x
            M2(i) = mean(B.^2);  % % M2 = <[x(t+tau)-x(t)]^2>|x(t)=x
            M4(i) = mean(B.^4);
        end
        A1(tau, :) = M1;
        A2(tau, :) = M2;
        A4(tau, :) = M4; % Note, to account for finite-tau corrections to diffusion consider A2(tau,:) = M2-M1.^2.
        ErrorM1(tau, :) = (M2 - M1.^2) ./ N;
        ErrorM2(tau, :) = (M4 - M2.^2) ./ N;
    end
else
    % Nadarya-Watson estimation of conditional moments
    
    h = 0.5 *(length(dat(1, :)))^(-1/4) * std(dat(1, :));
    K = @(x)(x>-sqrt(1)&x<sqrt(1)).*3/4 .*(3-x.^2); % Epanechnikov kernel
    %         h = 1.3 *(length(dat(1, :)))^( -1/ 4) * std(dat(1, :));
    %         K = @(x)(x>-sqrt(1) & x<sqrt(1)).*3 /4 .*(0.8-x.^2); % Epanechnikov kernel
    
    
    %         K=@(x)exp(-x.^2 ./2); % In case you wish to use Gaussian Kernel.
    for tau = Tau
        for i = 1:bins    %  Here, I use Matlab parallel computing toolbox. If you do not have it just use 'for loop' not 'parfor loop'.
            dat1 = dat;
            A = K((dat1(1, :) - C(i))./ h);  
            B = dat1(tau + 1, :) - dat1(1, :);  % x(t+tau)-x(t)
            ndx = ~isnan(B);
            A = A(ndx);
            B = B(ndx);
            N1 = sum(A .* B);
            N2 = sum(A .* B.^2);
            N4 = sum(A .* B.^4);   
            E = sum(A);
            M1(i) = N1 / E;  
            M2(i) = N2 / E;   
            M4(i) = N4 / E; 
        end
        A1(tau, :) = M1;
        A2(tau, :) = M2;
        A4(tau, :) = M4;
        ErrorM1(tau, :) = (M2 - M1.^2) ./ N;
        ErrorM2(tau, :) = (M4 - M2.^2) ./ N;
    end
end
WeightM1 = 1 ./ ErrorM1;
WeightM2 = 1 ./ ErrorM2;

% This is to make sure all weights are positive.
WeightM1(WeightM1 < 0 | isnan(WeightM1) | isinf(WeightM1)) = 0;
WeightM2(WeightM2 < 0 | isnan(WeightM2) | isinf(WeightM2)) = 0;

% Estimating the limit conditional moments/dt as dt-->0
D1 = zeros(1, bins);
D2 = D1;
D4 = D1;
Tau1 = dt .* Tau;
ErrorD1 = D1;
ErrorD2 = D1;
ErrorD4 = D1;

%preselect fit methods for efficiency 
%fit = fittype( 'poly1' );   
%opts = fitoptions( 'Method', 'LinearLeastSquares' );
% Calculation of Drift (D1), Diffusion (D2) and, D4
for i = 1:bins
    
    %     fit1=fitlm(Tau,A1(:,i),'weight',WeightM1(:,i));ErrorD1(i)=fit1.Coefficients.SE(2);
    %     fit2=fitlm(Tau,A2(:,i),'weight',WeightM2(:,i));ErrorD2(i)=fit2.Coefficients.SE(2);
    %     fit4=fitlm(Tau,A4(:,i),'weight',ones(1,length(Tau)));ErrorD4(i)=fit4.Coefficients.SE(2);
    %
    %     D1(i)=fit1.Coefficients.Estimate(2);
    %     D2(i)=1/2*fit2.Coefficients.Estimate(2);
    %     D4(i)=1/24*fit4.Coefficients.Estimate(2);
    
    % drfit
    fit1 = fitlm(Tau1, A1(:, i), 'weight', WeightM1(:, i));
    ErrorD1(i) = fit1.Coefficients.SE(2);
    D1(i) = fit1.Coefficients.Estimate(2); 
    
    % diffusion
    fit2 = fitlm(Tau1, A2(:, i) - (D1(i) .* Tau1').^2, 'weight', WeightM2(:, i));
    ErrorD2(i) = fit2.Coefficients.SE(2);
    D2(i) = 1 / 2 * fit2.Coefficients.Estimate(2); 
    
    % D(4)
    fit4 = fitlm(Tau1, A4(:, i), 'weight', ones(1, length(Tau1)));
    ErrorD4(i) = fit4.Coefficients.SE(2);
    D4(i) = 1 / 24 * fit4.Coefficients.Estimate(2);
end

if nargout == 1 
    D1 = struct('D1', D1, 'D2', D2, 'D4', D4, 'ErrorD1', ErrorD1, 'ErrorD2', ErrorD2, 'N', N, 'C', C, ...
        'options', struct('datasize', size(data), 'domain', [L, R], 'bins', bins, 'Tau', Tau, 'dt', dt, 'method', method, 'title', title));
end
end

function dat1 = makesegments(dat, len_segments)
%simple function to make segments of the time series dat
%we use this for bootstrapping
%
% example:
% makesegments(1:10,5)
%
% ans =
%
%      1     2     3     4     5     6
%      2     3     4     5     6     7
%      3     4     5     6     7     8
%      4     5     6     7     8     9
%      5     6     7     8     9    10
%      6     7     8     9    10   NaN
%      7     8     9    10   NaN   NaN
%      8     9    10   NaN   NaN   NaN
%      9    10   NaN   NaN   NaN   NaN
%     10   NaN   NaN   NaN   NaN   NaN
N = length(dat);
dat1 = nan(len_segments + 1, N);
dat1(1, :) = dat(1:N);
for i = 1:len_segments
    dat1(i + 1, 1:N - i) = dat(1 + i:end);
end
end

function [fitresult, gof] = LinearFit(a, b, w, ft, opts)
%CREATEFIT8(A,B,W)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : a
%      Y Output: b
%      Weights : w
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 05-Apr-2019 18:56:03


%% Fit: 'untitled fit 1'.
[xData, yData, weights] = prepareCurveData( a, b, w );

% Set up fittype and options.
if nargin < 4
    ft = fittype( 'poly1' );
end
if nargin < 5
    opts = fitoptions( 'Method', 'LinearLeastSquares' );
end
% ft = fittype( 'smoothingspline' );
% opts = fitoptions( 'Method', 'SmoothingSpline' );

opts.Weights = weights;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'b vs. a with w', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel a
% ylabel b
% grid on


end

function e = SlopeError(x, y, y_est, w)
n = length(x);
% w=w./sum(w);
e = sqrt(sum(w .* (y - y_est).^2) / ((n) * (sum(w .* (x - mean(x)).^2))));
end
