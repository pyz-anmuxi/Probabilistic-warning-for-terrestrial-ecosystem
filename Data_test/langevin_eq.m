classdef langevin_eq < handle
    %This class is used to analyze a reconstructed Langevin equation. It is
    %especially suitable for discretized D1 and D2 functions.
    %When you use smooth functions they are also descretized
    %Use this class to calculate/plot the propeties of the Langevin
    %equation:
    %Main methods:
    %   plot: used to plot the results
    %   pdf: calculated the stationary pdf
    %   mean_exit: calculates the mean exit time for all initial conditions
    %       of a domain
    %   survival: calculates the survival function for all initial
    %       conditions
    %   potential: the potential function based on D1
    %   potential_eff: the effective potential function, includes effects of
    %   D2
    %
    %Optionally, it uses chebfun: 
    %   Driscoll, T. A., N. Hale, and L. N. Trefethen. 2014. Chebfun Guide. Pafnuty Publication, Oxford.
    %
    properties
        D1 = []; %drift
        D2 = []; %diffusion
        DD2 = []; %derivative of diffusion
        equilibria = struct('x', {}, 'stable', false, 'BC', '', 'domain', []); %list of equilibria
        % x = state variable, stable = logical variable true if stable, BC boundary conditons RA AA AR or RR, 
        % domain - the basin of attraction of stable equilibrium, outer borders are -Inf Inf
        xtra %struct with extra information, for instance last results and reconstr information
        nx = 100; % default number of x for Fokker-Planck equations
        domain = [0 10]; % domain of the Langevin equation
        namex = 'x' % name of the x coordinate for plotting
        timeunit = ''; % time unit for plotting
    end
    
    methods
        function res = getinitfun(obj, atype, varargin)
            %  get an initial function for pdepe/bvp5v:
            %  usage:
            %  obj.getinitfun('normal_dom') - for runpde
            %  obj.getinitfun('normal',mean,std) - for runpde
            %  obj.getinitfun('dirichlet',x0) - for runpde
            %  obj.getinitfun('cosine',BC) - for mean_exit
          
            function u0 = dirichlet_delta_ic(x0, meshx) %vectorized
                dx = meshx(2) - meshx(1);
                [~, ix0] = min(abs(meshx - x0));
                u0 = zeros(size(meshx));
                ndx = abs(meshx - x0) < dx / 2 | (abs(meshx - x0) == dx / 2 & meshx > x0);
                u0(ndx) = 1 / dx;
                if ix0 == 1 || ix0 == length(meshx)
                    u0 = u0 * 2;
                end
            end
            %initial sine, make sure that the BC are correct
            %initfun = @(x)[cos(4 * x), -4 * sin(4 * x) ]   ;
            switch atype
                case 'normalize'
                    fun = varargin{1};
                    domain_ = obj.domain;
                    sumpdf = integral(fun, domain_(1), domain_(2));  
                    res = @(x)fun(x)/sumpdf;
                case 'normal_dom'
                    res = obj.getinitfun('normalize', @(x)exp(-(x - obj.domain(2) ./ 2).^2));
                case 'normal'
                    mean = varargin{1};
                    std = varargin{2};
                    res = @(x)pdf('normal',x,mean,std);
                case 'dirichlet'
                    x0 = varargin{1};
                    if length(varargin) > 1
                        meshx = varargin{2};
                    else
                        meshx = linspace(obj.domain(1), obj.domain(2), obj.nx);
                    end
                    u0 = dirichlet_delta_ic(x0, meshx);
                    interpolant = griddedInterpolant(meshx, u0);
                    res = @(x)interpolant(x);
                    % res = @(x)dirichlet_delta_ic(x, x0, meshx);
                case 'cosine;sine' %for mean_exit
                    % BC the boundary value code:
                    %  'AA' = left and right absorbing
                    %  'RA' = left reflecting and right absorbing
                    %  'AR' = left absorbing and right reflecting
                    domain_ = obj.domain;
                    BC = varargin{1};
                    starth = 0;
                    endh = 2 * pi;
                    if BC(1) == 'R'
                        %start at pi/2 to have correct BC
                        starth = starth + pi / 2;
                    end
                    if BC(2) == 'R'
                        %start at pi/2 to have correct BC
                        endh = endh + pi / 2;
                    end
                    per = (endh - starth) / (domain_(2) - domain_(1));
                    res = @(x)[cos(per * (x-domain_(1))+starth); -per * sin(per * (x-domain_(1))+starth) ];
            end
        end
    end
    
    methods (Access = protected)
        
% pdepe solve forward equation
        function result = solvepdf(obj, pdf0, varargin)
            %solve pdfs
            %used both by pdf and by runpdf
            if isempty(obj.D1) || isempty(obj.D2)
                error('langevin_eq:nomodel', 'No Langevin equation defined');
            end
            if nargin > 2 && ischar(pdf0) 
                varargin = [{pdf0} varargin];
                pdf0 = [];
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%% initial pdf0
            if nargin < 2 || isempty(pdf0)
                if isa(obj.D1, 'chebfun')
                    pdf0 = chebfun(@(x)exp(-(x - obj.domain(2) / 2)^2), obj.domain, 'vectorize');
                    pdf0 = pdf0 / sum(pdf0);
                else
                    pdf0 = obj.getinitfun('normal_dom');   % @(x)fun(x)/sumpdf    其中sumpdf=integral(fun, domain_(1), domain_(2))
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%% option
            if length(varargin) > 1
                options = struct(varargin{:});
            elseif length(varargin) == 1
                options = varargin{1};
            else
                options = [];
            end  
            if ~isfield(options, 'maxtime')
                options.maxtime = 1000;
            end
            if ~isfield(options, 'ntime')
                options.ntime = 1000;
            end
            if ~isfield(options, 'nx')
                options.nx = obj.nx;
            end

            L = min(obj.domain);
            R = max(obj.domain);
            m = 0;
            result.t = linspace(0, options.maxtime, options.ntime);
            x = linspace(L, R, options.nx);
               % Mesh = x;
            if isnumeric(pdf0) && length(pdf0) == 1  
                pdf0 = obj.getinitfun('dirichlet', pdf0, x);
                if isa(obj.D1, 'chebfun')
                    pdf0 = chebfun(pdf0, obj.domain, 'vectorize');
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%% pdepe function: fp_equat
            if isempty(obj.DD2)
                DD2_ = diff_fun(obj.D2); 
            else
                DD2_ = obj.DD2;
            end
            %we bind  the D1 D2 and DD2 functions to this anonymous function 
            if isa(obj.D1, 'chebfun')
                %using chebfuns is 100x slower, griddedInterpolant are really
                %efficient
                x1 = linspace(L, R, options.nx * 1000); %more points for precision
                d1 = griddedInterpolant(x1, obj.D1(x1));
                d2 = griddedInterpolant(x1, obj.D2(x1));
                dd2 = griddedInterpolant(x1, DD2_(x1));
                fp_equat = @(x,t,u,DuDx)pdex1pde(x, 0, u, DuDx ,d1,d2,dd2);   
            else
                fp_equat = @(x,t,u,DuDx)pdex1pde(x, 0, u, DuDx ,obj.D1, obj.D2, DD2_);  % 函数句柄
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%% option
            if ~isfield(options, 'RelTol') || isempty(options.RelTol)
                options.RelTol = 1E-3;
            end
            if ~isfield(options, 'AbsTol') || isempty(options.AbsTol)
                options.AbsTol = 1E-4;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%% boundary condition
            if ~isfield(options, 'BC') || isempty(options.BC)
                options.BC = 'RR';
            end
            switch options.BC
                case 'AA'
                    pdebc = @(~, ul, ~, ur, ~)deal(ul,0,ur,0);
                case 'AR'
                    pdebc = @(~, ul, ~, ~, ~)deal(ul,0,0,1);
                case 'RA'
                    pdebc = @(~, ~, ~, ur, ~)deal(0,1,ur,0);
                otherwise %default RR
                    pdebc = @(varargin)deal(0,1,0,1); %border conditions; pl, ql pr qr
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            sol = pdepe(m, fp_equat, pdf0, pdebc, x, result.t, options);
            % sol = pdepe(m,@pdex1pde,@pdex1ic,@pdex1bc,x,t);
            u = sol(:, :, 1);
            result.x = x;
            result.pdfs = u;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [c, f, s] = pdex1pde(x, ~, u, DuDx , D1, D2, DD2)
                %Forward Fokker Planck:
                %diff(P,t) = - diff (D1(x)*P(x,t),x) + diff(diff(D2(x)*P(x,t),x),x)
                %
                %general form for pdepe (see Matlab help)
                %c(x,t,u,dudx)*diff(u,t) = x^m * diff(x^m * f(x,t,u,dudx),x) + s(x,t,u,dudx)
                %
                %we need work out our double integral as that is not in pdepe (we define DD2= d(D2)/dx)
                %product rule: diff(D2(x)*P(x,t),x)= DD2 * P + D2 * dPdx
                %so our reworked equation is
                %diff(P,t) = diff (-D1 * P + DD2 * P + D2 * dPdx,x)
                %
                %this means we loose some terms:
                %m=0: c=1: s=0
                % keeping only:
                %diff(u,t) = diff(f,x)
                %where
                %f =-D1*P + DD2 * P + D2 * dPdx
                %   
                c = 1;
                f = -D1(x) .* u + DD2(x) .* u + D2(x) .* DuDx;
                s = 0;
            end
        end

        function [res, isconstant] = set_Dfun(obj, Dfun, avar)
            if nargin < 3
                avar = '';
            end
            isconstant = false;
            res = Dfun;
            if isa(Dfun, 'chebfun')
                obj.domain = Dfun.domain;
                if strcmp(avar, 'g')
                    res = 0.5 * res^2;
                end
            elseif ischar(Dfun)
                val = str2double(Dfun);
                if ~isnan(val)
                    [res, isconstant] = obj.set_Dfun(val, avar);
                    return;
                else
                    [res, isconstant] = obj.set_Dfun(evalin('base', sprintf('@(%s)%s', obj.namex, Dfun)), avar);
                end

            elseif isnumeric(Dfun) && numel(Dfun) == 1  
                isconstant = true;
                if strcmp(avar, 'g')
                    res = @(x)0.5*Dfun^2+zeros(size(x));
                else
                    res = @(x)Dfun+zeros(size(x));
                end

            elseif isa(Dfun, 'function_handle') 
                try
                    if strcmp(avar, 'g') 
                        Dfun = @(x)0.5*Dfun(x).^2;
                    end
                    %test whether the function is in vector notation or
                    %returning a constant.
                    s1 = Dfun(1:10);   
                    if numel(s1) == 1   
                        %constant
                        isconstant = true;
                        res = @(x)Dfun(x)+zeros(size(x));
                    elseif numel(s1) ~= 10 
                        %wrong size of resultr
                        res = @(x)simulate_vector(Dfun,x);
                    else
                        res = Dfun;
                    end
                catch err
                    %no vector notation supported, but simulate a
                    %vector
                    if strcmp(err.identifier, 'MATLAB:innerdim')
                        res = @(x)simulate_vector(Dfun,x);
                    end
                end
            end
        end
    end
 
    methods
        function obj = langevin_eq(varargin)
            %if chebfun is on the search path it will be used for the
            %functions (except if option 'nocheb' is used)
            %
            %obj = langevin_eq(res,'weightedspline') - where res is
            %the result of the function LangevinReconst, 'weightedspline'
            %fits a spline while weighing down points that are uncertain.

            %is the default if the curvefitting toolbox is present (else
            %spline)
            %obj = langevin_eq(res,'spline') - uses standard spline
            %interpolation (not recommended)
            %obj = langevin_eq(res,'weightedspline','nocheb') - do not use
            %chebfun functions
            %
            %obj = langevin_eq('D1',@func1,'D2',@funct2,'domain',[min
            %max]) - function handle mode of this class, give
            %function handles for D1 and D2
            %D2 can also be a constant (just enter a numeric value)
            %
            if nargin > 0
                ndx = cellfun(@ischar, varargin);  
                options = varargin(ndx); 
                
%%%%%%%%%%%%%%%%%%%%%%%%%%
                if nargin > 0 && isfield(varargin{1}, 'ErrorD2') 
                    %result of LangevinReconst splining the D1 and D2 functions determined with the Langevin Approach
                    %(aka Langevin Reconstruction)
                    results = varargin{1};
                    if ~isfield(results, 'options')
                        results.options = options;
                    end
                    if isfield(results, 'bootstrap')
                        obj.xtra.bootstrap = results.bootstrap;
                        results = rmfield(results, 'bootstrap'); 
                    end
                    
                    L1 = results.C(1);
                    R1 = results.C(end);
                    
                    %reconstructed Langevin equation
                    if ((nargin == 1 && exist('fittype', 'file')) || any(strcmp(options, 'weightedspline')))  % exist: 检测im中的变量是否存在  any: 检测矩阵中是否有非零元素
                        %weighted spline fitting  Weighing down of
                        %uncertain points
                        if ~exist('fittype', 'file')
                            error('langevin_eq:toolboxmissing', 'The curvefitting toolbox is required for the ''weightedspline'' option')
                        end
                        ft = fittype( 'smoothingspline' );
                        %
                   %%%%%%%%%%%%%fit D1
                        if all(results.ErrorD1 == 0) 
                            w1 = ones(size(results.ErrorD1));
                        else
                            w1 = 1 ./ results.ErrorD1;
                        end
                        [xData, yData, weights1] = prepareCurveData(results.C, results.D1, w1 );  
                        opts = fitoptions( 'Method', 'SmoothingSpline', 'SmoothingParam', 0.9982070110205727, 'Weights', weights1);  % fitoptions: 
                        [d1, results.gofD1] = fit( xData, yData, ft, opts );
                        %
                    %%%%%%%%%%%%% fit D2
                        if all(results.ErrorD2 == 0)
                            w2 = ones(size(results.ErrorD1));
                        else
                            w2 = 1 ./ results.ErrorD2;
                        end
                        w2(results.D2 - results.ErrorD2 <= 0) = 0; %I assign 0 weight to the very few bins for D2 when D2<0
                        [xData, yData, weights1] = prepareCurveData(results.C, results.D2, w2 );
                        opts = fitoptions( 'Method', 'SmoothingSpline', 'SmoothingParam', 0.9982070110205727, 'Weights', weights1);
                        [d2, results.gofD2] = fit( xData, yData, ft, opts );
                        
                        
                        if exist('chebfun', 'file') && ~any(strcmp(options, 'nocheb'))
                            l = linspace(L1, R1, 2^13);
                            D1 = chebfun(d1(l), [L1 R1], 'equi', 5000);  
                            D2 = chebfun(d2(l), [L1 R1], 'equi', 5000);
                            obj.set('D1', D1, 'D2', D2);  % Inf×1 chebfun
                        else
                            obj.set('D1', @(x)reshape(d1(x),size(x)), 'D2', @(x)reshape(d2(x),size(x)), 'domain', [L1 R1]);  % B = reshape(A,sz) 使用大小向量sz重构A以定义size(B)
                        end
                        
                    elseif nargin == 1 || any(strcmp(options, 'spline')) %default if curvefitting toolbox is missing
                        d1 = griddedInterpolant(results.C, results.D1, 'spline'); 
                        d2 = griddedInterpolant(results.C, results.D2, 'spline');
                        if exist('chebfun', 'file') && ~any(strcmp(options, 'nocheb'))
                            D1 = chebfun(@(x)d1(x), [L1 R1]);
                            D2 = chebfun(@(x)d2(x), [L1 R1]);
                            obj.set('D1', D1, 'D2', D2);
                        else
                            obj.set('D1', @(x)d1(x), 'D2', @(x)d2(x), 'domain', [L1 R1]);
                        end
                    end
                    obj.xtra.reconstr = results;
                % 2 default arguments D1 and D2  
                elseif nargin == 2 && ~ischar(varargin{1})
                    obj.set('D1', varargin{1}, 'D2', varargin{2});
                % openrealmod 
                else
                    obj.set(varargin{:}); 
                end
                % eguilibrium
                if isempty(obj.equilibria)
                    obj.equilibria = obj.find_equilibria('deterministic'); %can be changed to 'effective' if you use ueff
                end
                
            end
        end
 
        function res = simple_stoch_exit(obj, data, varargin)
            %one-step (or n-step) conditional distribution 
            %default nstep=1
            %
            
            function res = SojournTimes(data, dt, stable_equilibrium, width)
                repelor = stable_equilibrium.domain(~isinf(stable_equilibrium.domain));
                if numel(repelor) > 1
                    error('not yet implemented');
                end
                res.periods_in=cell(size(data,2),1);
                for m = 1:size(data, 2)
                    s=[];
                    data1 = data(:, m);
                    if stable_equilibrium.x < repelor
                        T = 0;
                        %discard the first part of the time series if we start below the
                        %stable equilibrium
                        Idx = find(data1 > stable_equilibrium.x - width, 1);
                        if Idx > 1
                            Idx = Idx - 1;
                        end
                        data1 = data1(Idx:end);
                        while ~isempty(T)
                            Idx = find(data1 < stable_equilibrium.x + width, 1);
                            %Idx=find(data>e-width & data<e+width,1);
                            data1 = data1(Idx:end);
                            T = find(data1 > repelor, 1);
                            if ~isempty(T)
                                s = [s T];
                                data1 = data1(T:end);
                            end
                        end
                    else
                        T = 0;
                        %discard the first part of the time series if we start above the
                        %stable equilibrium
                        Idx = find(data1 < stable_equilibrium.x + width, 1);
                        if Idx > 1
                            Idx = Idx - 1;
                        end
                        data1 = data1(Idx:end);
                        while ~isempty(T)
                            Idx = find(data1 > stable_equilibrium.x - width, 1);
                            %Idx=find(data>e-width & data<e+width,1);
                            data1 = data1(Idx:end);
                            T = find(data1 < repelor, 1);
                            if ~isempty(T)
                                s = [s T];
                                data1 = data1(T:end);
                            end
                        end
                    end
                    res.periods_in{m} = s .* dt;
                end 
                if size(data, 2)==1
                    res.periods_in=res.periods_in{m};
                end
                res.equilibrium = stable_equilibrium;
            end
            
            %             function res = direct_exit(datat, datax, domain)
            %                 %simple way of determining exit time from data or simulations:
            %                 % we first (optionally) apply a moving average to smooth the data set. 
            %                 % we then determine the mean exit times based on transitions past a certain
            %                 % threshold
            %                 res.periods_in = [];
            %                 res.domain = domain;
            %                 for m = 1:size(datax, 2)
            %                     higheq = datax(:, m) > domain(1) & datax(:, m) < domain(2);
            %                     shiftsup = datat(diff(higheq) == 1);
            %                     shiftsdown = datat(diff(higheq) == -1);
            %                     period_up = nan(size(shiftsdown));
            %                     for i = 1:length(shiftsup)
            %                         k = find(shiftsdown > shiftsup(i), 1);
            %                         %we get nan if none is found
            %                         if ~isempty(k)
            %                             period_up(i) = shiftsdown(k) - shiftsup(i);
            %                         end
            %                     end
            %                     res.periods_in = [res.periods_in; period_up];
            %                 end
            %                 ndx = ~isnan(res.periods_in);
            %                 res.periods_in = res.periods_in(ndx);
            %                 % res.period_down = period_down;
            % 
            %             end

            args = struct(varargin{:});
            if ~isfield(args, 'windowsize')
                args.windowsize = 1;
            end
            if ~isfield(args, 'dt')
                args.dt = obj.xtra.reconstr.options.dt;
            end
            if ~isfield(args, 'nreplicates')
                args.nreplicates = 100;
            end
            if ~isfield(args, 'width')
                args.width = 0;
            end
            if ~isfield(args, 'MaxStep')
                args.MaxStep = 0.001;
            end
            if ~isfield(args, 'NonNegative')
                args.NonNegative = [];
            end
            if args.windowsize <= 1
                trend = data(:);
            else
                trend = movmean(data(:), args.windowsize);
            end
            % datat = linspace(0, numel(data) * args.dt, numel(data))';
            direct_exits = cell(1, obj.nbasin);
            % simul_exits = cell(1, obj.nbasin);
            stableeq = find([obj.equilibria.stable]);
            res.equilibria = [obj.equilibria(stableeq).x];
            for i1 = 1:length(stableeq)
                direct_exits{i1} = SojournTimes(data, args.dt, obj.equilibria(stableeq(i1)), args.width);
                %direct_exits{i1} = direct_exit(datat(:), trend, obj.equilibria(stableeq(i1)).domain);
                %     simul_exits{i1} = direct_exit(simuldatat, simultrend, obj.equilibria(stableeq(i1)).domain);
            end
            res.dataexits = direct_exits;
            %res.simulexits = simul_exits;
            
            
            if ~isfield(args, 'simuldata')
                disp('Simulating');
                sizedata = numel(data);
                [simuldata, fcn] = obj.simulate( linspace(0, (sizedata - 1) * args.dt, sizedata)', obj.draw_pdf(args.nreplicates, 1), struct('MaxStep', args.MaxStep, 'NonNegative', args.NonNegative));
                res.simuldata = fcn;
                simuldatat = simuldata.t;
                simuldata = simuldata.y;
                if args.windowsize <= 1
                    simultrend = simuldata;
                else
                    simultrend = movmean(simuldata, args.windowsize);
                end
            end

            res.options = args;
            datat = linspace(0, numel(data) * args.dt, numel(data))';
            direct_exits = cell(1, obj.nbasin);
            simul_exits = cell(1, obj.nbasin);
            stableeq = find([obj.equilibria.stable]);
            res.equilibria = [obj.equilibria(stableeq).x];
            for i1 = 1:length(stableeq)
                direct_exits{i1} = SojournTimes(trend, args.dt, obj.equilibria(stableeq(i1)), args.width);
                %direct_exits{i1} = direct_exit(datat(:), trend, obj.equilibria(stableeq(i1)).domain);
                simul_exits{i1} = SojournTimes(simultrend, args.dt, obj.equilibria(stableeq(i1)), args.width);
                %simul_exits{i1} = direct_exit(simuldatat, simultrend, obj.equilibria(stableeq(i1)).domain);
            end
            res.dataexits = direct_exits;
            res.simulexits = simul_exits;
            %             figure
            %             subplot(2, 1, 1)
            %             stairs(sort(res.simulexits{1}.periods_in), cumsum(sort(res.simulexits{1}.periods_in)) / sum(res.simulexits{1}.periods_in))
            %             hold on
            %             stairs(sort(res.dataexits{1}.periods_in), cumsum(sort(res.dataexits{1}.periods_in)) / sum(res.dataexits{1}.periods_in))
            %             subplot(2, 1, 2)
            %             stairs(sort(res.simulexits{2}.periods_in), cumsum(sort(res.simulexits{2}.periods_in)) / sum(res.simulexits{2}.periods_in))
            %             hold on
            %             stairs(sort(res.dataexits{2}.periods_in), cumsum(sort(res.dataexits{2}.periods_in)) / sum(res.dataexits{2}.periods_in))
        end
   
        function res = condition_step_distr(obj, data, varargin)
            %one-step (or n-step) conditional distribution 
            %default nstep=1
            %
            function [distr, sizes] = getdistrib(thedata, stepdata, perc)
                if nargin < 3
                    perc = [];
                end
                distr = cell(args.bins, 1);
                sizes = zeros(args.bins, 1);
                bincenters = linspace(args.domain(1), args.domain(2), args.bins);
                mesh = linspace(args.domain(1), args.domain(2), args.nmesh);
                binsiz = (bincenters(2) - bincenters(1)) / 2;
                for l = 1:args.bins
                    F = nan(size(thedata, 2), args.nmesh);
                    siz = nan(size(thedata, 2), 1);
                    for j = 1:size(thedata, 2)
                        ndx1 = thedata(:, j) > bincenters(l) - binsiz & thedata(:, j) < bincenters(l) + binsiz;
                        data1 = stepdata(ndx1, j);
                        siz(j) = length(data1);
                        if ~isempty(data1)
                            [F1, ~] = ksdensity(data1(:), mesh);
                            F(j, :) = F1 ./ (trapz(mesh, F1));
                        else
                            F(j, :) = 0;
                        end
                    end
                    if ~isempty(perc)
                        distr{l} = quantile(F, perc);
                        sizes(l) = median(siz);
                    else
                        sizes(l) = siz;
                        distr{l} = F;
                    end
                end
            end
            args = struct(varargin{:});
            if ~isfield(args, 'bins')
                args.bins = obj.xtra.reconstr.options.bins;
            end
            if ~isfield(args, 'nmesh')
                args.nmesh = obj.nx;
            end
            if ~isfield(args, 'dt')
                args.dt = obj.xtra.reconstr.options.dt;
            end
            if ~isfield(args, 'nreplicates')
                args.nreplicates = 100;
            end
            if ~isfield(args, 'MaxStep')
                args.MaxStep = 0.001;
            end
            if ~isfield(args, 'NonNegative')
                args.NonNegative = [];
            end
            if ~isfield(args, 'nstep')
                args.nstep = 1;
            end
            if ~isfield(args, 'domain')
                args.domain = obj.domain;
            end
            if ~isfield(args, 'percentiles')
                args.percentiles = [0.275 0.5 0.975];
            end
            data = data(:);
            sizedata = size(data, 1);
            data1 = data(args.nstep + 1:end);
            stepdata = data(1:end - args.nstep); %second row is the time lag
            res.bincenters = linspace(args.domain(1), args.domain(2), args.bins);
            res.mesh = linspace(args.domain(1), args.domain(2), args.nmesh);
            [res.datadistr, res.sizes] = getdistrib(data1, stepdata);
            disp('Simulating');
            [simuldata, fcn] = obj.simulate( linspace(0, (sizedata - 1) * args.dt, sizedata)', ...
                obj.draw_pdf(args.nreplicates, 1), struct('MaxStep', args.MaxStep, 'NonNegative', args.NonNegative));
            res.simuldata = fcn;
            simuldata1 = simuldata.y(args.nstep + 1:end, :);
            stepsimuldata = simuldata.y(1:end - args.nstep, :);
            res.simuldatadistr = getdistrib(simuldata1, stepsimuldata, args.percentiles);
            n1 = 10;
            if n1 > args.nreplicates
                n1 = args.nreplicates;
            end
            res.exampledata = simuldata.y(:, 1:n1); %at most ten example data sets 
            res.options = args;
            sumdiff = 0;
            sumsiz = sum(res.sizes);
            for l = 1:length(res.sizes)
                simuldist = res.simuldatadistr{l};
                simuldist = simuldist(2, :);
                datadist = res.datadistr{l};
                sumdiff = sumdiff + res.sizes(l) / sumsiz * sum(abs(simuldist - datadist));
            end
            res.sumdiff = sumdiff;
        end
        
% Calculate the equilibrium point, stability, and boundary conditions、domain
        function equilibria1 = find_equilibria(obj, mode)
            %find the equilibria and basin of attractions
            % can be improved as it finds sometimes bubbles
            %if mode is 'deterministic'
            % the equilibria are just roots of D1
            %if mode is 'effective' (default)
            % the effect of D2 is accounted
            if nargin == 1
                mode = 'effective';
            end
            if strncmpi(mode, 'd', 1) 
                %naively D1
                d1 = obj.D1;
                mode = 'deterministic';
            else
                %-derivative of the effective potential (chain rule on log)
                if isa(obj.D1, 'chebfun')
                    d1 = obj.D1 ./ obj.D2 - diff(obj.D2) ./ obj.D2;
                else
                    dd2 = diff_fun(obj.D2);
                    d1 = @(x)obj.D1(x) ./ obj.D2(x) - dd2(x) ./ obj.D2(x);
                end
                mode = 'effective';
            end
            
            Eqs = roots_fun(d1, obj.domain);   % Solving for the equilibrium point
            DD1 = diff_fun(d1);   
            isstable = DD1(Eqs) < 0; 
            
            Eqms = struct('x', {}, 'stable', [], 'BC', [], 'domain', []);  

            for i = 1:length(Eqs)
                Eqms(i) = struct('x', Eqs(i), 'stable', isstable(i), 'BC', '', 'domain', []);
            end

            for i = 1:length(Eqms)
                if isstable(i)
                    BC = 'RR';
                    dom = [-Inf Inf];
                    % [Eqms(2).x Eqms(3).x]
                    if i > 1
                        BC(1) = 'A';
                        dom(1) = Eqms(i - 1).x;
                    end
                    % [Eqms(1).x Eqms(2).x]
                    if i < length(Eqms)
                        BC(2) = 'A';
                        dom(2) = Eqms(i + 1).x;
                    end
                    Eqms(i).domain = dom;
                    Eqms(i).BC = BC;
                end
            end
            
            if nargout == 0   
                fprintf('Equilibria (%s):\n', mode)
                for i = 1:length(Eqms)
                    if Eqms(i).stable
                        stab = '  stable';
                    else
                        stab = 'unstable';
                    end
                    fprintf('%s equilibrium: %g   %s  %s\n', stab, Eqms(i).x, Eqms(i).BC, mat2str(Eqms(i).domain));
                end
            else
                equilibria1 = Eqms;
            end
        end
        
% ode function
        function res = odefun(obj, dt)
            %creates an ode function to be used in Euler integration
            %it can also 
            L = obj.domain(1);
            R = obj.domain(2);
            Mesh = linspace(L, R, 1000);
            %griddedInterpolant is really fast
            d1 = griddedInterpolant(Mesh, obj.D1(Mesh), 'linear');  
            d2 = griddedInterpolant(Mesh, obj.D2(Mesh), 'linear');
            %the noise part  is multiplied with sqrt(t) and divided by dt, 
            %in Euler it is multiplied again with dt, so it is a correct Euler Maruyama 
            %
            %if x<L, x=L
            %if x>R, x=R
            %D1(x) + sqrt(2*dt*D2(x))*randn(1)/dt
            res = @(t,y)d1(min(R,max(L,y)))+sqrt(2 .* dt .* d2(min(R,max(L,y)))) .* randn(size(y)) /dt;
        end
        
% Solve the SDE to generate time series data.
        function [res, fcn] = simulate(obj, tspan, y0, options)
            %simulate - time integration of the Langevin equation using
            % Euler-Maruyama
            % 
            % Arguments:
            % tspan = time span of the integration or all t values
            % y0    = the initial condition (default 0.01). y0 can be a vector, the 
            %         model will then be integrated for each y0 (efficient
            %         due to vector notation)
            % options = a structure with options for the integration, used
            %         fields: MaxStep (or StepSize) - the step size
            %                 NonNegative - force non-negative see Matlab
            % Outputs:
            % res   struct with the data sets (res.t= time; res.y = data)
            % fcn   function handle that can reproduce the data (use for
            %       efficient storage)
            % 
            %     
            if nargin < 2
                tspan = 1:1000;
            end
            if nargin < 3
                y0 = 0.01;
            end
            if nargin < 4
                options = odeset([], 'AbsTol', []);
            end
            if isfield(options, 'MaxStep') && ~isempty(options.MaxStep)
                dt = options.MaxStep;
            else
                if isfield(options, 'StepSize') && ~isempty(options.StepSize)
                    dt = options.StepSize;
                else
                    dt = 0.01;
                end
            end
            %the noise part is divided by dt, in Euler it is multiplied
            %again, to it is a correct Euler Mayaruma
            %save random seed
            r1 = rng;  
            %some prework
            odefun = obj.odefun(dt);
            options = odeset(options); 
            options.MaxStep = dt;
            fcn = @()simulate_rng(r1,odefun,tspan,y0,options);
            %[res.t, res.y] = euler(odefun, tspan, y0, options);
            res = fcn();
        end
      
% 
        function res = draw_pdf(obj, siz1, siz2)
            %draw_pdf - draw random number from the stationary pdf
            %use this to get a stabilized initial condition
            %
            % res=obj.draw_pdf(siz1,siz2)
            if ~isfield(obj.xtra, 'results') || ~isfield(obj.xtra.results, 'pdf')
                thepdf = obj.pdf; 
            else
                thepdf = obj.xtra.results.pdf;
            end
            %make cumulative density function
            cdf = cumtrapz(thepdf.x, thepdf.PDF);  
            if nargin < 2
                siz1 = 1;
            end
            if nargin < 3
                siz2 = siz1;
            end
            %interpolate from the inverted cumulative density function
            res = interp1(cdf, thepdf.x, rand(siz1, siz2));  
        end
        
        function adjusttimescale(obj, ts)
            %optionnally you can adjust the time scale of the Langevin
            %equation by multiplying D1 and D2 with a factor
            if ~isfield(obj.xtra.reconstr, 'ts')
                obj.D1 = obj.D1 * ts;
                obj.D2 = obj.D2 * ts;
                obj.xtra.reconstr.ts = ts;
                if isfield(obj.xtra.reconstr, 'D1')
                    obj.xtra.reconstr.D1 = obj.xtra.reconstr.D1 * ts;
                    obj.xtra.reconstr.D2 = obj.xtra.reconstr.D2 * ts;
                    obj.xtra.reconstr.ErrorD1 = obj.xtra.reconstr.ErrorD1 * ts;
                    obj.xtra.reconstr.ErrorD2 = obj.xtra.reconstr.ErrorD2 * ts;
                end
            else
                error('Timescale was already changed');
            end
        end

        function result = potential(obj)
            %calculate the potential function
            if isempty(obj.D1) || isempty(obj.D2)
                error('langevin_eq:nomodel', 'No Langevin equation defined');
            end
            if isa(obj.D1, 'chebfun')
                result.u = -cumsum(obj.D1); % Clearly, here we do not need to care about time scale change.
            else
                h = integralx_fun(obj.D1, obj.domain(1), obj.domain(2));
                result.u = @(x)-h(x);
            end
            result.dom = linspace(obj.domain(1), obj.domain(2), obj.nx);
        end

        function result = potential_eff(obj)
            %calculate the effective potential. This includes the effects
            %of diffusion on the stationary pdf
            if isempty(obj.D1) || isempty(obj.D2)
                error('langevin_eq:nomodel', 'No Langevin equation defined');
            end
            if isa(obj.D1, 'chebfun')
                result.ueff = cumsum(-obj.D1 / obj.D2) + log(obj.D2); % Clearly, here we do not need to care about time scale change.
            else
                h1 = @(x)-obj.D1(x)/obj.D2(x);
                h = integralx_fun(h1, obj.domain(1), obj.domain(2));
                result.ueff = @(x)h(x)+log(obj.D2(x));
            end
            result.dom = linspace(obj.domain(1), obj.domain(2), obj.nx);
        end
        

        function res = prob_exit(obj, domain_, BC, options)
            % prob_exit (numerically) determines the probability to exit
            % left or right
            %
            % obj.mean_exit(ibasin, PDF, options)
            %
            % Inputs are:
            %
            % domain= range of conditions
            % BC the boundary value code:
            %  'AA' = left and right absorbing
            %  'RA' = left reflecting and right absorbing
            %  'AR' = left absorbing and right reflecting
            %
            % options (optional) the boundary value problem options (see
            % bvpopt) + additional fields:
            %    nx = number of points in x direction
            %    initfunc = function handle describing the initial guess of
            %    the boundary value problem
            %  can also be entered as field, value pairs (case sensitive)
            %
            % Outputs is a structure with fields:
            % domain: the domain
            % BC: boundary condition
            % P_left:  the probability of exit from left boundary.
            % P_right: the probability of exit from right boundary.
            %
            function res = getinitfun(domain_, BC)
                %initial sine, make sure that the BC are correct
                %initfun = @(x)[cos(4 * x), -4 * sin(4 * x) ]   ;
                starth = 0;
                endh = 2 * pi;
                if BC(1) == 'R'
                    %start at pi/2 to have correct BC
                    starth = starth + pi / 2;
                end
                if BC(2) == 'R'
                    %start at pi/2 to have correct BC
                    endh = endh + pi / 2;
                end
                per = (endh - starth) / (domain_(2) - domain_(1));
                res = @(x)[cos(per * (x-domain_(1))+starth); -per * sin(per * (x-domain_(1))+starth) ];
            end
            res.domain = domain_;
            res.BC = BC;
            if ~strcmp(BC, 'AA')
                if BC(1) == 'R'
                    res.P_left = 0;
                else
                    res.P_left = 1;
                end
                if BC(2) == 'R'
                    res.P_left = 0;
                else
                    res.P_left = 1;
                end
            else

                if isa(obj.D1, 'chebfun')
                    d1 = obj.D1;
                    d1.domain = domain_;
                    d2 = obj.D2;
                    d2.domain = domain_;
                    LP = chebop(domain_); % LP is operator for the probability of exit
                    LP.op = @(x, P) d1 .* diff(P) + d2 .* diff(P, 2);
                    LP.lbc = @(P)P - 1;
                    LP.rbc = @(P)P;
                    P = LP \ (0);
                    res.P_left = P;

                    LP = chebop(domain_); % LP is operator for the probability of exit
                    LP.op = @(x, P) d1 .* diff(P) + d2 .* diff(P, 2);
                    LP.lbc = @(P)P;
                    LP.rbc = @(P)P - 1;
                    P = LP \ (0);
 
                    res.P_right = P;
                else
                    % function handle mode
                    %
                    %Backward Fokker Planck equation:
                    %D1(x0)*T'+D2(x0 )*T''=-1 
                    %
                    %define:
                    %T2=T';
                    %D1(x0)*T2+D2(x0)*T2'=-1
                    %
                    %yields the following 2D ODE system
                    %T2'=-(D1(x0)*T2 + 1)/D2(x0)
                    %T'=T2;
                    if nargin < 4
                        options = [];
                    end
                    if ~isfield(options, 'nx')
                        options.nx = obj.nx;
                    end
                    if ~isfield(options, 'RelTol') || isempty(options.RelTol)
                        options.RelTol = 1E-7;
                    end
                    if ~isfield(options, 'AbsTol') || isempty(options.RelTol)
                        options.RelTol = 1E-6;
                    end
                    if ~isfield(options, 'initfun') || isempty(options.initfun)
                        options.initfun = {getinitfun(domain_, 'AR'), getinitfun(domain_, 'RA')};
                    elseif ~iscell(options.initfun)
                        options.initfun = {options.initfun, options.initfun};
                    end
                    xpoints = linspace(res.domain(1), res.domain(2), options.nx);
                    %We determine the probabilities of leaving the
                    %left border with het BVP equation:
                    %0 = D1(x0) * P' + D2(x0) * P'';
                    %
                    %define:
                    %P2  =P';
                    %D1(x0) *P2 + D2(x0) * P2' = 0
                    %
                    %yields the following 2D ODE system
                    %P'=P2;
                    %P2'=-(D1(x0)*P2)/D2(x0)
                    %
                    odefunP = @(x,P)[P(2,:);-obj.D1(x) .* P(2,:) ./ obj.D2(x) ];

                    %left BC P = 1, right BC P = 0
                    BCfunleft = @(ya,yb)[ya(1)-1; yb(1)];
                    solinit = bvpinit(xpoints, options.initfun{1});
                    %determine the left probabilities
                    sol = bvp4c(odefunP, BCfunleft, solinit, options);
                    res.P_left = deval(sol, xpoints); %only if we want to interpolate between the mesh
                    res.P_left = res.P_left(1, :);
                    BCfunright = @(ya,yb)[ya(1); yb(1)-1];
                    solinit = bvpinit(xpoints, options.initfun{2});
                    sol = bvp4c(odefunP, BCfunright, solinit);
                    res.P_right = deval(sol, xpoints); %only if we want to interpolate between the mesh
                    res.P_right = res.P_right(1, :);
                end
            end
        end
        
% bvp4 calculate exit_time (edge ??value problem)
        function res = mean_exit(obj, ibasin, PDF, varargin)
            % mean_exit (numerically) solves the backward Fokker-Planck equation
            % and then calculates many quantities.
            %
            % obj.mean_exit(ibasin, PDF, options)
            %
            % Inputs are:
            %
            % ibasin = the number of the basin or 'all' or a struct with
            % the fields "domain" and "BC"
            % Where domain is the domain and BC the boundary value code:
            %  'AA' = left and right absorbing
            %  'RA' = left reflecting and right absorbing
            %  'AR' = left absorbing and right reflecting
            %
            % PDF (optional) is a struct describing the stationary probability density
            % function or empty, see output of the method "pdf"
            %
            % options (optional) the boundary value problem options (see
            % bvpopt) + additional fields:
            %    nx = number of points in x direction
            %    initfunc = function handle describing the initial guess of
            %    the boundary value problem
            %  can also be entered as field, value pairs (case sensitive)
            %
            % Outputs is a structure with fields:
            %
            % T: is a (chebfun) function which tells us the average time to escape the
            % domain (see exit curves in Figures 2,3,5 in the main text).
            % WT: is the average time to escape the domain. It implements the calculations of Formula 4 in the main text.
            % P: calculates the probability of exit from left or right boundaries. Note that, it only makes sense if
            % the boundary conditions are of AAleft or AAright types (we did not use it in this paper).
            % 
            % use obj.plot(res) to plot the results
            %
  
            if nargin < 2
                ibasin = 'all';
            end
            if nargin < 3
                PDF = [];
            end
            if nargin > 4
                options = struct(varargin{:});
            elseif nargin == 4
                options = varargin{1};
            else
                options = bvpset; 
            end
            
           
            if ischar(ibasin) && strcmp(ibasin, 'all')  
                res = cell(1, obj.nbasin);   % obj.nbasin
                for i = 1:length(res)
                    res{i} = obj.mean_exit(i, PDF);  % PDF = []
                end
                obj.xtra.results.mean_exit = res;
                return;
            end
          
            if isstruct(ibasin) && length(ibasin) > 1
                %if ibasin is a struct with equilibria remove basins without BC
                ndx = ~cellfun('isempty', {ibasin(:).BC});
                ibasin = ibasin(ndx);
                res = cell(1, length(ibasin));
                for i = 1:length(res)
                    res{i} = obj.mean_exit(ibasin(i), PDF);
                end
                obj.xtra.results.mean_exit = res;
                return;
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% （ 'RA'，‘AR’），domain
            if isstruct(ibasin)
                domain2 = ibasin.domain;
                domain2(domain2 < obj.domain(1)) = obj.domain(1);
                domain2(domain2 > obj.domain(2)) = obj.domain(2);
                BC = ibasin.BC;
            else
                nbasin1 = obj.nbasin;
                if ibasin > nbasin1
                    error('langevin_eq:mean_exit', 'basin number too large');
                end
                stableeq = find([obj.equilibria.stable]); 
                eqnr = stableeq(ibasin);
                domain2 = obj.equilibria(eqnr).domain;
                domain2(domain2 < obj.domain(1)) = obj.domain(1);
                domain2(domain2 > obj.domain(2)) = obj.domain(2);
                BC = obj.equilibria(eqnr).BC;
            end

            if ~any(strcmp(BC, {'RA', 'AR', 'AA'}))
                error('langevin_eq:mean_exit', 'boundary condition "%s" unknown', BC);
            end
            
            domain_ = domain2;
            res.domain = domain_;
            res.BC = BC;
            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% initfun
            if ~isfield(options, 'initfun') || isempty(options.initfun)
                %initial sine, make sure that the BC are correct
                %initfun = @(x)[cos(4 * x), -4 * sin(4 * x) ]   ;
                starth = 0;
                endh = 2 * pi;
                if BC(1) == 'R'
                    %start at pi/2 to have correct BC
                    starth = starth + pi / 2;
                end
                if BC(2) == 'R'
                    %start at pi/2 to have correct BC
                    endh = endh + pi / 2;
                end
                per = (endh - starth) / (res.domain(2) - res.domain(1));
                %                 if isa(obj.D1,'chebfun')
                %                     x=chebfun('x');
                %                     initfunT=[cos(per * (x-res.domain(1))+starth); -per * sin(per * (x-res.domain(1))+starth) ] ;
                %                 else    
                initfunT = @(x)[cos(per * (x - res.domain(1)) + starth); -per * sin(per * (x - res.domain(1)) + starth) ]   ;
                %               end
                %do not add initfun to options, as we use them also for
                %prob_exit
            else
                initfunT = options.initfun;
            end       
%%%%%%%%%%%%%%%%%%%%%% 
            if isa(obj.D1, 'chebfun')
                if ~isempty(PDF)
                    if ~isa(PDF, 'chebfun')
                        PDF = chebfun(PDF.PDF', obj.domain, 'equi', 5000);
                    end
                    PDF.domain = domain2;
                end
                d1 = obj.D1;
                d1.domain = domain_;
                d2 = obj.D2;
                d2.domain = domain_;
                %setting up chebop boundary value problem

                LL = chebop(domain_);
                LL.op = @(x, T) d1 .* diff(T) +d2 .* diff(T, 2);
                if BC(1) == 'R'
                    LL.lbc = @(T)diff(T);
                else
                    LL.lbc = @(T)T;
                end
                if BC(2) == 'R'
                    LL.rbc = @(T)diff(T);
                else
                    LL.rbc = @(T)T;
                end
                %solving T
                res.T = solvebvp(LL , -1);
                T2 = chebfun(res.T, domain2);
                
                if ~isempty(PDF)
                    res.WT = sum(PDF * T2) / sum(PDF);
                else
                    res.WT = [];
                end
                res.P = [];

                if strcmp(BC, 'AA')
                    res1 = obj.prob_exit(domain_, BC);
                    LL = chebop(domain_);
                    LL.op = @(x, T) d1 .* diff(T) + d2 .* diff(T, 2);
                    LL.lbc = @(T)T;
                    LL.rbc = @(T)T;
                    % T = LL \ (-res1.P_left);
                    [T, sol] = solvebvp(LL , -res1.P_left);
                    if sol.error > 0.001
                        xpoints = linspace(domain_(1), domain_(2), 200);
                        P_l = griddedInterpolant(xpoints, res1.P_left(xpoints));
                        odefunT = @(x,T)[T(2)+zeros(size(x));  -(obj.D1(x) .* T(2) + P_l(x)) ./ obj.D2(x)];
                        BCfun = @(ya,yb)[ya(1) yb(1)];
                        one = chebfun(1, domain_);
                        v0 = [one 0 * one];
                        T1 = bvp4c(odefunT, BCfun, v0);
                        res.T_left = chebfun(T1(:, 1) ./ res1.P_left, 'splitting', 'on');
                    else
                        res.T_left = chebfun(T ./ res1.P_left, 'splitting', 'on');
                    end

                    T2 = chebfun(res.T_left, domain2, 'splitting', 'on');
                    if ~isempty(PDF)
                        res.WT_left = sum(PDF * T2) / sum(PDF);
                    else
                        res.WT_left = [];
                    end
                    res.P_left = res1.P_left;

                    LL = chebop(domain_);
                    LL.op = @(x, T) d1 .* diff(T) + d2 .* diff(T, 2);
                    LL.lbc = @(T)T;
                    LL.rbc = @(T)T;
                    %T = LL \ (- res1.P_right);
                    [T, sol] = solvebvp(LL , -res1.P_right);
                    if sol.error > 0.001
                        xpoints = linspace(domain_(1), domain_(2), 200);
                        P_l = griddedInterpolant(xpoints, res1.P_right(xpoints));
                        odefunT = @(x,T)[T(2)+zeros(size(x));  -(obj.D1(x) .* T(2) + P_l(x)) ./ obj.D2(x)];
                        BCfun = @(ya,yb)[ya(1) yb(1)];
                        one = chebfun(1, domain_);
                        v0 = [one 0 * one];
                        T1 = bvp4c(odefunT, BCfun, v0);
                        res.T_right = chebfun(T1(:, 1) ./ res1.P_right, 'splitting', 'on');
                    else
                        res.T_right = chebfun(T ./ res1.P_right, 'splitting', 'on');
                    end
                    %plot(T);

                    T2 = chebfun(res.T_right, domain2, 'splitting', 'on');
                    if ~isempty(PDF)
                        res.WT_right = [];
                    else
                        res.WT_right = sum(PDF * T2) / sum(PDF);
                    end
                    res.P_right = res1.P_right;
                end
 %%%%%%%%%%%%%%%%%%%%% 
            else
                % function handle mode
                %
                %Backward Fokker Planck equation:
                %D1(x0)*T'+D2(x0 )*T''=-1 
                %
                %define:
                %T2=T';
                %D1(x0)*T2+D2(x0)*T2'=-1
                %
                %yields the following 2D ODE system
                %T2'=-(D1(x0)*T2 + 1)/D2(x0)
                %T'=T2;
                if ~isfield(options, 'nx')
                    options.nx = obj.nx;
                end
                if ~isfield(options, 'RelTol') || isempty(options.RelTol)
                    options.RelTol = 1E-7;
                end
                if ~isfield(options, 'AbsTol') || isempty(options.RelTol)
                    options.RelTol = 1E-6;
                end
                options.Vectorized = 'on';
                xpoints = linspace(res.domain(1), res.domain(2), options.nx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BCfun
                %Reflecting means that T'==0
                %Absorbing means that T==0
                %this BCfun function should return the RESIDUAL
                %ya is left yb is right boundary
                switch BC
                    case 'RA'
                        BCfun = @(ya,yb)[ya(2) yb(1)];
                    case 'AR'
                        BCfun = @(ya,yb)[ya(1) yb(2)];
                    case 'AA'
                        BCfun = @(ya,yb)[ya(1) yb(1)];
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  solinit          
                solinit = bvpinit(xpoints, initfunT); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% The BVP solver returns the structure 'sol'. 
                %for exit time T(x)
                try  % 执行语句并捕获产生的错误
                    odefunT = @(x,y)[y(2,:);  -(obj.D1(x) .* y(2,:) + 1) ./ obj.D2(x)];
                    sol = bvp4c(odefunT, BCfun, solinit, options);
                    Tx = deval(sol, xpoints);  
                catch err
                    disp(err.message); 
                    Tx = NaN;
                end
                res.T = Tx(1, :);
                if ~isempty(PDF)
                    pdf = griddedInterpolant(PDF.x, PDF.PDF);
                    PDF = pdf(xpoints);
                    res.WT = sum(PDF .* res.T) ./ sum(PDF);
                else
                    res.WT = [];
                end
                
                if strcmp(BC, 'AA')
                    %first we need to determine the probabilities of leaving the
                    %left border with BVP equation:
                    res1 = obj.prob_exit(domain_, BC, options);
                    %mean exit through left:
                    res.P_left = res1.P_left;
                    P_l = griddedInterpolant(xpoints, res1.P_left);
                    odefunT = @(x,T)[T(2,:)+zeros(size(x));  -(obj.D1(x) .* T(2,:) + P_l(x)) ./ obj.D2(x)];

                    sol = bvp4c(odefunT, BCfun, solinit, options);
                    res.T_left = deval(sol, xpoints);
                    res.T_left = res.T_left(1, :) ./ res.P_left;

                    res.P_right = res1.P_right;
                    P_r = griddedInterpolant(xpoints, res.P_right);
                    odefunT = @(x,T)[T(2,:);  -(obj.D1(x) .* T(2,:) + P_r(x)) ./ obj.D2(x)];
                    sol = bvp4c(odefunT, BCfun, solinit, options);
                    res.T_right = deval(sol, xpoints);
                    res.T_right = res.T_right(1, :) ./ res.P_right;
                end
            end
            obj.xtra.results.mean_exit = res;
        end

% Numerical solution of the Fokker-Planck equations
        function result = runpdf(obj, varargin)
            % This code (numerically) solves the Fokker-Planck equation.
            % 
            % obj.solvepdf(PDF0,options)
            % Inputs are:
            % PDF0: is the initial distribution of states (is a chebfun function or function handle 
            % or a scalar number (Dirichlet delta function is then assumed).
            %
            % 
            % options (optional) the pdepe options (see pdepe, for instance RelTol and AbsTol)
            %  + additional fields:
            %    maxtime = duration of the simulation
            %    ntime = number of time outputs (minimal=3) (default 3)
            %    nx = number of points in x direction (default obj.nx)
            %
            %  can also be entered as field, value pairs (case sensitive)
            %
            % Output is a structure with fields:
            %   x = mesh of x values
            %   t = time output values
            %   pdfs = pdfs in time
            % 
            result = solvepdf(obj, varargin{:});
            obj.xtra.results.pdfs = result;
        end

% Two methods for solving steady-state probability distribution
        function result = pdf(obj, pdf0, varargin)
% method 1：This code (numerically) solves the Fokker-Planck equation.
            % 
            % obj.pdf(PDF0,options)
            %
            % Inputs are:
            % PDF0: is the initial distribution of states (is a chebfun function or function handle 
            % or a scalar number (Dirichlet delta function is then assumed).
            %
            % options (optional) the pdepe options (see pdepe, for instance RelTol and AbsTol)
            %  + additional fields:
            %    maxtime = duration of the simulation
            %    ntime = number of time outputs (minimal=3) (default 3)
            %    nx = number of points in x direction (default obj.nx)
            %
            %  can also be entered as field, value pairs (case sensitive)
            %
            % Output is a structure with fields:
            %   x = mesh of x values
            %   t = time output values
            %   pdfs = all pdfs in time
            %   PDF = last (stationary) pdf
 
            %initial 
            if nargin < 2
                pdf0 = [];
            end
            %option
            if nargin > 3
                options = struct(varargin{:});
            elseif nargin == 3
                options = varargin{1};  
            else
                options = [];
            end
            if ~isfield(options, 'maxtime') 
                options.maxtime = 4000;
            end
            if ~isfield(options, 'ntime')
                options.ntime = 3;
            end
            % 
            result = obj.solvepdf(pdf0, options);  % 
            result.PDF = result.pdfs(end, :);  % 
            
% method 2：another way of determining the PDF using the effective potential
            %integral D1/D2:
            %
            %
            if isa(obj.D1, 'chebfun')
                try  
                    h2 = 1 / obj.D2 * exp(cumsum(obj.D1 / obj.D2));
                    
                    %normalize to surface = 1
                    xpoints = linspace(obj.domain(1), obj.domain(2), 4000);
                    intg = trapz(xpoints, h2(xpoints)); 
                    %normalized pdf
                    result.PDF_eq = 1 / obj.D2 * exp(cumsum(obj.D1 / obj.D2)) / intg;
                catch
                    result.PDF_eq = [];
                end
            else
                h1 = @(x)obj.D1(x)./obj.D2(x);
                h = integralx_fun(h1, obj.domain(1), obj.domain(2));   
                %
                %1/D2 * exp(int(D1/D2,x))
                h2 = @(x)1 ./obj.D2(x).*exp(h(x));

                %normalize to surface = 1
                xpoints = linspace(obj.domain(1), obj.domain(2), 4000);
                intg = trapz(xpoints, h2(xpoints)); 
                %normalized pdf
                result.PDF_eq = @(x)1 ./obj.D2(x).*exp(h(x))./intg;
            end
            obj.xtra.results.pdf = result;
        end
        

        function res = survival(obj, ibasin, varargin)
            % This code (numerically) calculates the survival function (Formula 1 in the main text) of the 'Langevin system' for all initial states.
            % You can find full mathematical details in Appendix S3 which is from the following book
            % Schuss, Z. Theory and applications of stochastic processes: an analytical approach. Vol. 170 (Springer Science & Business Media, 2009).
            %
            % Inputs are:
            %
            % ibasin = the number of the basin or 'all' or a struct with
            % the fields "domain" and "BC"
            % Where domain is the domain and BC the boundary value code:
            %  'AA' = left and right absorbing
            %  'RA' = left reflecting and right absorbing
            %  'AR' = left absorbing and right reflecting
            %
            % options (optional) the pdepe options (RelTol and AbsTol for instance)
            %    + additional fields:
            %    time_range = range of time 
            %    ntime = number of points for output of time axis
            %    nx = number of points in x direction
            %    initfunc = function handle describing the initial
            %    conditions (default = all 1)
            % options can also be entered as field, value pairs (case
            % sensitive)
            % for instance:
            % obj.survival('all','ntime',100,'AbsTol',1E-6)
            % 
            % IMPORTANT NOTE: 
            % For a much better accuracy we recommend NOT to use a regular mesh for mesh_time. Instead, we
            % recommend to use a chebychev mesh point which has 2^n points. For instance, 
            % to generate a chebychev mesh point on the interval [0 300] of length 2^16
            % simply type mesh_time = chebpts([0 300],2^15);
            %
            % Output is:
            % structure with fields:
            %    t = output points in time
            %    x = meshpoints of x
            %    survival = output of pdepe
            %  you can use obj.plot(res) to plot the output structure
            %
            if nargin < 2
                ibasin = 'all';
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%% options：time_range，ntime，nx
            if nargin > 3
                options = struct(varargin{:});
            elseif nargin == 3
                options = varargin{1};  % [0 400]
            else
                options = [];
            end
            if isnumeric(options) && ~isempty(options) %alternatively you can only give the time_range
                options = struct('time_range', options);
            end
            if ~isfield(options, 'time_range')
                options.time_range = []; % empty is automatic (takes an extra run however)
            end
            if ~isfield(options, 'ntime') 
                options.ntime = 500;
            end
            if ~isfield(options, 'nx')
                options.nx = obj.nx;
            end
 
            
            if ischar(ibasin) && strcmp(ibasin, 'all')
                res = cell(1, obj.nbasin);  % obj.nbasin: 
                for i = 1:length(res)  
                    res{i} = obj.survival(i, options); 
                end
                obj.xtra.results.survival = res;
                return;
            end
            
            if isstruct(ibasin)
                eqn = ibasin;
            else
                nbasin1 = obj.nbasin;  % obj.nbasin
                if ibasin > nbasin1
                    error('langevin_eq:survival', 'basin number too large');
                end
                stableeq = find([obj.equilibria.stable]); 
                eqn = obj.equilibria(stableeq(ibasin)); 
            end
            
            % 保证区间在obj.domain内
            domain_ = eqn.domain;
            domain_(domain_ < obj.domain(1)) = obj.domain(1);  % left < domain时
            domain_(domain_ > obj.domain(2)) = obj.domain(2);  % right > domain时
            
            % 
            switch eqn.BC  % 
                    %absorbing = ul,0 (or ur,0)  reflecting = 0,1
                case 'AA'
                    BCfun = @(~, ul, ~, ur, ~)deal(ul, 0 , ur, 0); %"deal" deals the values to the 4 requested outputs
                case 'RA'
                    BCfun = @(~, ~, ~, ur, ~)deal(0, 1, ur, 0);
                case 'AR'
                    BCfun = @(~, ul, ~, ~, ~)deal(ul, 0, 0, 1);
                otherwise
                    error('langevin_eq:bc', 'boundary conditions undefined');
            end
            
            if ~isfield(options, 'Stats')
                options.Stats = 'on';
            end
            if ~isfield(options, 'RelTol') || isempty(options.RelTol)
                options.RelTol = 1E-6;
            end
            if ~isfield(options, 'AbsTol') || isempty(options.AbsTol)
                options.AbsTol = 1E-7;
            end
            if ~isfield(options, 'initfun') || isempty(options.initfun)
                options.initfun = @(x)1; 
            end
            
            if isempty(obj.DD2)
                DD2_ = diff_fun(obj.D2);
            else
                DD2_ = obj.DD2;
            end
            
            %
            %we bind  the functions D1 D2 and DD2_ in the anonymous(匿名) function
            if isa(obj.D1, 'chebfun')
                %else this is much slower
                x = linspace(min(domain_), max(domain_), 100 * options.nx);
                d1 = griddedInterpolant(x, obj.D1(x));
                d2 = griddedInterpolant(x, obj.D2(x));
                dd2 = griddedInterpolant(x, DD2_(x));
                surv_equat = @(x,t,u,DuDx)survpde(x, 0, u, DuDx, d1, d2, dd2);
            else
                surv_equat = @(x,t,u,DuDx)survpde(x, 0, u, DuDx ,obj.D1, obj.D2, DD2_);
            end
            
            m = 0;
            if isempty(options.time_range)
                options.time_range = [0 500];
                %                 %do a test run to find a reasonable range for the time
                %                 %(without having too many zeros - does
                %                 %not work very well
                options1 = options;
                %                 options1.time_range=[0 1000];
                options1.ntime = 500;
                testrun = obj.survival(1,options1);
                meansurv = max(testrun.survival,[],2);
                options.time_range = [0 find(meansurv<0.1,1)+1];
            end
            
            % 
            if length(options.time_range) == 2
                res.t = linspace(options.time_range(1), options.time_range(2), options.ntime)';  % ntime = 500 
            elseif length(options.time_range) == 1
                res.t = linspace(0, options.time_range(1), options.ntime)';
            else
                res.t = linspace(0, options.time_range(end), options.ntime)';
            end
            %res.t = options.time_range; % t = linspace(0,T,N_t);  % IMPORTANT: here t is not the usual time we use for IVPs (initial value problems).
            %Our problem here is IBVP (initial boundary value problem). Therefore considering only 3 mesh points for t does not work.
            
            % 
            res.x = linspace(min(domain_), max(domain_), options.nx);  % nx = 100 

            %options = odeset('Stats', 'on', 'RelTol', 1e-6, 'AbsTol', 1e-7);
            res.survival = pdepe(m, surv_equat, options.initfun, BCfun, res.x, res.t, options);
            obj.xtra.results.survival = res;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            function [c, f, s] = survpde(x, ~, ~, DuDx, D1, D2, DD2)
                %Survival PDE:
                %diff(S,t) = D1(x) diff(S(x,t),x) + D2(x)*diff(diff(S(x,t),x),x)
                %
                %general form for pdepe (see Matlab help)
                %>>> c(x,t,u,dudx)*diff(u,t) = x^m * diff(x^m * f(x,t,u,dudx),x) + s(x,t,u,dudx)
                %simplify:
                %m=0;c=1
                %                
                %assume f = D2*dSdx
                %remember D2 can be dependent on x (multiplicative noise)!
                %
                %so product rule: diff(D2 * dSdx,x)=D2 * diff(dSdx,x)+DD2 * dSdx
                %we only need D2*diff(dSdx,x) so we just subtract DD2*dSdx from the s term
                %m=0:c=1:
                %diff(u,t) = diff(f(x,t,u,dudx),x) + s(x,t,u,dudx)
                %f = D2 * dSdx
                %s = D1 * dSdx - DD2 * dSdx
                %
                c = 1;
                f = D2(x) * DuDx;
                s = (D1(x) - DD2(x)) * DuDx;
            end

        end
        
% The number of attraction basins, i.e., the number of stable equilibrium points.
        function n = nbasin(obj)
            %number of basins of attraction
            n = sum([obj.equilibria.stable]);
        end
        

        function set(obj, varargin)
            if nargin > 1
                if ~isstruct(varargin{1})
                    args = struct(varargin{:});   
                else
                    args = varargin{1};
                end
            end

            if isfield(args, 'domain')
                obj.domain = args.domain;
                if isa(obj.D1, 'chebfun')
                    obj.D1.domain = obj.domain;
                end
                if isa(obj.D2, 'chebfun')
                    obj.D2.domain = obj.domain;
                end
            end
            
            if isfield(args, 'timeunit')
                obj.timeunit = args.timeunit;
            end
            if isfield(args, 'namex')
                obj.namex = args.namex;
            end
            if isfield(args, 'nx')
                obj.nx = args.nx;
            end
            if isfield(args, 'f')
                obj.D1 = obj.set_Dfun(args.f);
            end
            if isfield(args, 'g')
                obj.D2 = obj.set_Dfun(args.g, 'g');
            end

            if isfield(args, 'D1')
                obj.D1 = obj.set_Dfun(args.D1);
            end
            if isfield(args, 'D2')
                [obj.D2, isconstant] = obj.set_Dfun(args.D2);
                if isconstant
                    % the derivative is zero
                    obj.DD2 = obj.set_Dfun(0);
                end
            end
            if isfield(args, 'DD2')
                obj.DD2 = obj.set_Dfun(args.DD2);
            end
            if isfield(args, 'equilibria')
                obj.equilibria = args.equilibria;
            end
            if isfield(args, 'nx')
                obj.nx = args.nx;
            end
        end
        
% plot exit_time\survice\pst\Ueff
        function [hax, res] = plot(obj, varargin)
            %plot the results of the analyses, two ways of running:
            % 1) use the output of one of the methods as input of plot:
            %     >> obj.plot(obj.pdf) - to plot the stationary pdf
            % 2) use text: (if necessary the analyses is done with default
            %     parameters)
            %     for instance:
            %     >> figure;
            %     >> obj.plot('D1');
            %
            % parameters:
            %     'basinprob'     plot the probabilities of being in each
            %                     attraction basin vs. time
            %     'D1'            plot the drift function
            %     'D2'            plot the diffusion function
            %     'DD2'           plot the derivative of the diffusion function
            %     'equilibria'    plot the equilibria in the D1 plot
            %     'exit_distrib'  plot the probability density function of
            %                     exit times for a certain initial state
            %     'exitprob'      probability of exit from each boundary in a
            %                     AA boundary condition
            %     'mean_exit'     mean exit time of the basin
            %     'pdf'           plot the stationary pdf
            %     'potential'     potential
            %     'potential_eff' effective potential
            %     'survival'      plot the survival and median exit time
            %     'survival_func' plot the survival function for an inital
            %                     state
            %     'vertical_basin_boundaries'  indicate the basin boundaries 
            %                     in the current figure


            function res = getbootresults(bootstrap)
                C = bootstrap.results(1).reconstr.C;
                Ts = zeros(numel(bootstrap.results), numel(C));
                WTs = zeros(numel(bootstrap.results), numel(bootstrap.results(1).results.mean_exit));
                for i1 = 1:numel(bootstrap.results)
                    T = zeros(size(C));
                    for j1 = 1:numel(bootstrap.results(i1).results.mean_exit)
                        me = bootstrap(1).results(i1).results.mean_exit{j1};
                        ndx = C >= me.domain(1) & C <= me.domain(2);
                        WTs(i1, j1) = me.WT;
                        if isa(me.T, 'chebfun')
                            T(ndx) = me.T(C(ndx));
                        else
                            T(ndx) = me.T(ndx);
                        end
                    end
                    Ts(i1, :) = T;
                end
                if ~isfield(bootstrap.options, 'bias')
                    bias = NaN;
                else
                    bias = bootstrap.options.bias;
                end
                res = struct('C', C, 'Ts', Ts, 'WTs', WTs, 'CL_WT', quantile(WTs, [0.025, 0.975]), 'CL_Ts', quantile(Ts, [0.025, 0.975]), 'bias', bias);
            end
% plot mean_exit
            function hax1 = plot_mean_exit(T, domain, WT, acolor, nx)
                if isa(T, 'chebfun')
                    l1 = linspace(min(domain), max(domain), nx);
                    hax1 = plot(l1, T(l1), '-', 'LineWidth', 1.8, 'Color',[137/256 23/256 23/256]);
                else
                    l1 = linspace(min(domain), max(domain), numel(T));  
                    hax1 = plot(l1, T, 'LineWidth', 1.8, 'Color',[137/256 23/256 23/256]);    % 绘制退出时间（绿色的线）
                end
                hold on
                if ~isempty(WT) 
                    plot([l1(1), l1(end)], [WT, WT], '--', 'LineWidth', 2.5, 'Color', [156/256 204/256 101/256]);
                end
            end
            
            if isempty(obj.D1) || isempty(obj.D2)
                error('langevin_eq:nomodel', 'No Langevin equation defined');
            end
            
            hax1 = [];
            res = [];
            
            if nargin == 1
                if isfield(obj.xtra, 'results')
                    f = fieldnames(obj.xtra.results);
                    for i = 1:length(f)
                        res = obj.xtra.results.(f{i});
                        if isfield(res, 'figno')
                            figure(res.figno)
                        else
                            figure;
                        end
                        obj.plot(res);
                    end
                else
                    hax1 = obj.plot('D1');
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif iscell(varargin{1})  
                for i = 1:length(varargin{1})
                    hax1(i) = obj.plot(varargin{1}{i}, varargin{2:end});  
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            elseif ischar(varargin{1}) 
                args = varargin(2:end);  
                switch varargin{1}
                    case {'D1', 'D2' , 'DD2', 'D4'}
                        %extra options:
                        %timescale: simply scales yaxis
                        %npoints: number of interpolation points
                        %ylabel: specify the y label
                        %noboot: true/false suppress plotting of bootstrap
                        %error: true/false suppress plotting of errors
                        args = struct(args{:});
                        if strcmp(varargin{1}, 'D1')
                            x = obj.D1;
                        elseif strcmp(varargin{1}, 'D2')
                            x = obj.D2;
                        elseif strcmp(varargin{1}, 'DD2')
                            if isempty(obj.DD2)
                                x = diff_fun(obj.D2);
                            else
                                x = obj.DD2;
                            end
                        elseif strcmp(varargin{1}, 'D4')
                            x = griddedInterpolant(obj.xtra.reconstr.C, obj.xtra.reconstr.D4, 'spline');
                        end
                        %                         ploterror=true;
                        %                         ts=1;
                        %                         args.verticalbar=true;
                        if ~isfield(args, 'timescale')
                            args.timescale = 1;
                        end
                        if ~isfield(args, 'noboot')
                            args.noboot = false;
                        end
                        if ~isfield(args, 'npoints')
                            args.npoints = 50;
                        end
                        if ~isfield(args, 'ylabel')
                            if strcmp(varargin{1}, 'D1')
                                args.ylabel = 'Drift D1';
                            elseif strcmp(varargin{1}, 'DD2')
                                args.ylabel = 'd(D2)/dx';
                            else
                                args.ylabel = 'Diffusion D2';
                            end
                        end

                        dom = linspace(obj.domain(1), obj.domain(2), args.npoints);
 
                        if isfield(obj.xtra, 'bootstrap') && ~args.noboot
                            results = obj.xtra.reconstr;
                            hold on
                            %vectorbased fast:
                            recs = [obj.xtra.bootstrap(1).results(:).reconstr];
                            boots = vertcat(recs.(varargin{1}));
                            range = quantile(boots, [0.025, 0.975]);
                            fill([results.C , fliplr(results.C)], args.timescale * [range(1, :), fliplr(range(2, :))], [130/255 177/255 255/255], 'EdgeColor', 'none');
                        elseif ~(~isfield(obj.xtra, 'reconstr') || isempty(obj.xtra.reconstr)) && (~isfield(args, 'error') || args.error)
                            hold on
                            results = obj.xtra.reconstr;
                            if strcmp(varargin{1}, 'D1')
                                %hax1 = errorbar(results.C, results.D1, results.ErrorD1, '.b');
                                fill([results.C , fliplr(results.C)], [args.timescale * results.D1 - results.ErrorD1, fliplr(args.timescale * results.D1 + results.ErrorD1)], [130/255 177/255 255/255], 'EdgeColor', 'none')
                            elseif strcmp(varargin{1}, 'D2')
                                fill([results.C , fliplr(results.C)], [args.timescale * results.D2 - results.ErrorD2, fliplr(args.timescale * results.D2 + results.ErrorD2)], [130/255 177/255 255/255], 'EdgeColor', 'none')
                            end
                        end
                        
                        hax1 = plot(dom, args.timescale.* x(dom), 'color', [30/255 72/255 186/255], 'linewidth', 1.5, 'tag', varargin{1});
                        xlim(obj.domain);
                        xlabel(obj.namex);
                        ylabel(args.ylabel)
                        
                        if strcmp(varargin{1}, 'D1')
                            yline(0, 'k-');
                        end
                        if ~isfield(args, 'verticalbar') || args.verticalbar
                            obj.plot('vertical_basin_boundaries');
                        end
                        %                     case 'D4'
                        %                         %plot the Relative error of the Langevin
                        %                         %reconstruction
                        %                         if ~isfield(args, 'timescale')
                        %                             args.timescale = 1;
                        %                         end
                        %                         if isfield(obj.xtra, 'bootstrap')
                        %                             results = obj.xtra.reconstr;
                        %                             hold on
                        %                             %boots = zeros(numel(obj.xtra.bootstrap(1).results), numel(results.C));
                        %                             recs = [obj.xtra.bootstrap(1).results(:).reconstr];
                        %                             boots = vertcat(recs.D4);
                        %                             %                             for i = 1:numel(obj.xtra.bootstrap(1).results)
                        %                             %                                 boots(i, :) = obj.xtra.bootstrap(1).results(i).reconstr.D4;
                        %                             %                             end
                        %                             range = quantile(boots, [0.025, 0.975]);
                        %                             fill([results.C , fliplr(results.C)], [range(1, :), fliplr(range(2, :))], [1 0.8 0.8], 'EdgeColor', 'none');
                        %                         end
                        %                         hax1 = plot(obj.xtra.reconstr.C, args.timescale .* obj.xtra.reconstr.D4, '-r');
                        %                         hold on
                        %                         ylabel('D_{4}')
                        %                         xlabel(obj.namex);
                    case 'vertical_basin_boundaries'
                        hold on
                        unstableeqs = find(~[obj.equilibria.stable]);
                        for i = 1:length(unstableeqs)
                            hax1 = xline(obj.equilibria(unstableeqs(i)).x, 'k--');
                            %plot(obj.equilibria(unstableeqs(i)).x+zeros(1,2),ylim1,'k--')
                        end
                    case 'potential_eff'
                        hax1 = obj.plot(obj.potential_eff, args{:});
                    case 'potential'
                        hax1 = obj.plot(obj.potential, args{:});
                    case 'pdf'
                        if ~isfield(obj.xtra, 'results') || ~isfield(obj.xtra.results, 'pdf')
                            obj.pdf;
                        end
                        hax1 = obj.plot(obj.xtra.results.pdf);
                    case 'survival'
                        hold on
                        if nargin > 2 || ~isfield(obj.xtra, 'results') || ~isfield(obj.xtra.results, 'survival')
                            obj.survival('all', args{:});
                        end
                        hax1 = obj.plot(obj.xtra.results.survival);
                    case 'basinprob'
                        if nargin > 2
                            result = args{1};
                        else
                            if nargin == 2 && ~isfield(obj.xtra, 'results') || ~isfield(obj.xtra.results, 'pdfs')
                                obj.pdfs;
                            end
                            result = obj.xtra.results.pdfs;
                        end
                        nt = size(result.pdfs, 1);
                        result.t = result.t(1:nt);
                        eqs = obj.equilibria([obj.equilibria(:).stable]);
                        for i = 1:length(eqs)
                            eqs(i).domain = min(result.x(end), max(result.x(1), eqs(i).domain));
                        end
                        Ps = zeros(nt, length(eqs));

                        for j = 1:length(eqs)
                            ran = result.x < max(eqs(j).domain) & result.x > min(eqs(j).domain);
                            for i = 1:nt
                                Ps(i, j) = trapz(result.x(ran), result.pdfs(i, ran));
                            end
                        end
                        if ~isempty(Ps)
                            hax1 = plot(result.t, Ps);
                        end
                        xlabel('Time (t)')
                        ylabel('Probability')
                        leg = cell(length(eqs), 1);
                        for j = 1:length(eqs)
                            leg{j} = sprintf('%s in [%g,%g]', obj.namex, min(eqs(j).domain), max(eqs(j).domain));
                        end
                        legend(leg);
                        ylim([0, 1]);
                        xlim([result.t(1), result.t(end)]);
                    case 'exitprob' %args1 is mean_exit results of AA problem
                        if nargin > 2
                            result = args{1};
                        else
                            if nargin == 2 && ~isfield(obj.xtra, 'results') || ~isfield(obj.xtra.results, 'mean_exit')
                                obj.mean_exit;
                            end
                            result = obj.xtra.results.mean_exit;
                        end
                        if isfield(result, 'T') || isfield(result, 'T_left') %mean exit time
                            xpoints = linspace(result.domain(1), result.domain(2), length(result.P_left));
                            if isfield(result, 'T_left')
                                if isa(result.P_left, 'chebfun')
                                    plot(xpoints, result.P_left(xpoints), 'b');
                                    hold on
                                    plot(xpoints, result.P_right(xpoints), 'm');
                                else
                                    plot(xpoints, result.P_left, 'b')
                                    hold on
                                    plot(xpoints, result.P_right, 'm');
                                end
                            else
                                plot(xpoints, ones(size(xpoints)), 'b');
                            end
                            ylim([0 1])
                            xlabel(obj.namex)
                            ylabel('Probability of exit P')
                        end
                    case {'survival_func', 'exit_distrib'} %args1 is x0 args2 = survival results
                        hold on
                        x0 = args{1};
                        if (length(args) < 2) || ~(isstruct(args{2}) || iscell(args{2}))
                            if isfield(obj.xtra, 'results') && isfield(obj.xtra.results, 'survival')
                                results = obj.xtra.results.survival;
                            else
                                results = obj.survival('all', args{2:end});
                            end
                        else
                            results = args{2};  % cell(1*3)---struct  
                        end
                        
                        if iscell(results)
                            for i = 1:length(results)
                                [hax, res1] = obj.plot(varargin{1}, x0, results{i});
                                if ~isempty(res1)
                                    res = res1;
                                end
                            end
                        else
                            if x0 > min(results.x) && x0 < max(results.x) 
                                n = find(results.x >= x0, 1);  
                                surv = results.survival(:, n(1));   
                                res.surv = surv;
                                res.t = results.t;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                if strcmp(varargin{1}, 'survival_func')
                                    hax1 =fill([results.t(:); results.t(end); results.t(1)], [surv; 0; 0], [135 185 83]/256);  % [135 185 83]/256   [255 236 76]
                                    ylim([0 1]);
                                    ylabel('Survival S of high state')
                                end
                                
                                if strcmp(varargin{1}, 'exit_distrib')
                                    if isa(obj.D1, 'chebfun')
                                        f = chebfun(surv, [results.t(1) results.t(end)], 'equi', 10000);
                                        f = -diff(f);
                                    else
                                        f = @(x)interp1(results.t,surv,x,'spline');
                                        f1 = diff_fun(f);
                                        f = @(x)-f1(x);
                                    end
                                    res.density = f(results.t);
                                    hax1 = fill([results.t(:); results.t(end); results.t(1)], [res.density(:); 0; 0], [135 185 83]/256); % [135 185 83]/256   [255 236 76]/256
                                    ylabel('Density of exit times')
                                    ylim([0 Inf])
                                end
                                xlim([min(results.t), max(results.t)])
                                xlabel(sprintf('Time %s', obj.timeunit))
                            end
                        end
                    case 'mean_exit'
                        hold on
                        if isempty(args)
                            args = {'all', obj.pdf};
                        end
                        if isfield(args{1}, 'PDF')
                            args = [{'all'}, args];
                        end
                        if nargin > 2 || ~isfield(obj.xtra, 'results') || ~isfield(obj.xtra.results, 'mean_exit')
                            obj.mean_exit(args{:});
                        end
                        if isfield(obj.xtra, 'bootstrap')
                            %add the bootstrapped percentiles
                            hold on
                            res = getbootresults(obj.xtra.bootstrap(1));
                            fill([res.C , fliplr(res.C)], [res.CL_Ts(1, :), fliplr(res.CL_Ts(2, :))], [0.75 1 1], 'EdgeColor', 'none');
                            % range = quantile(res.WTs(:, [1 2]), [0.025 0.25 0.5 0.75 0.975]);
                            titl = [];
                            s = [];
                            disp(obj.xtra.bootstrap(1).bootmethod);
                            for i = 1:size(res.WTs, 2)
                                titl = sprintf('%sCL Basin %d\trelrange %d\t', titl, i, i);
                                CL = quantile(res.WTs(:, i), [0.025 0.975]);
                                med = median(res.WTs(:, i));
                                s = sprintf('%s[%g %g]\t%g\t', s, CL(1), CL(2), (CL(2) - CL(1)) / med);
                            end
                            fprintf('%s\n', titl);
                            fprintf('%s\n', s);
                        end
                        hax1 = obj.plot(obj.xtra.results.mean_exit);
                    case 'equilibria'
                        stab = [obj.equilibria([obj.equilibria.stable]).x];
                        unstab = [obj.equilibria(~[obj.equilibria.stable]).x];
                        plot(stab, zeros(size(stab)), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', 'k');
                        hold on
                        hax1 = plot(unstab, zeros(size(unstab)), 'ko', 'MarkerSize', 8, 'MarkerFaceColor', [1 1 1]);
                end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            else % varargin{1}-struct，survival、exit time：100*500数据
                
                result = varargin{1};
                args = struct(varargin{2:end});
                

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % mean exit time             
                if isfield(result, 'T') || isfield(result, 'T_left') 
                    if isfield(result, 'T_left')
                        plot_mean_exit(result.T_left, result.domain, [], 'b', obj.nx);
                        hold on;
                        plot_mean_exit(result.T_right, result.domain, [], 'm', obj.nx);
                    end
                    hax1 = plot_mean_exit(result.T, result.domain, result.WT, [137/256 23/256 23/256], obj.nx);  % 调用plot_mean_exit
                    xlim(obj.domain);
                    yl = get(gca, 'ylim');
                    if yl(1) < 0
                        set(gca, 'ylim', [0 yl(2)]);
                    end
                    xlabel(strcat('Initial state',32, obj.namex));  % 空格码：32
                    if isempty(obj.timeunit)
                        ylabel('Mean exit time T')
                    else
                        ylabel(sprintf('Mean exit time T (%s)', obj.timeunit))
                    end
                    if ~isfield(args, 'verticalbar') || args.verticalbar
                        obj.plot('vertical_basin_boundaries');
                    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % effective potential
                elseif isfield(result, 'ueff') 
                    U = result.ueff(result.dom);
                    U = U - min(U) + 0.1 * (max(U) - min(U));
                    hax1 = fill([result.dom, result.dom(end), result.dom(1)], [U, 0, 0], [130/255 177/255 255/255]);
                    xlim(obj.domain);
                    xlabel(obj.namex);
                    ylabel('Effective potential')
                    if ~isfield(args, 'verticalbar') || args.verticalbar
                        obj.plot('vertical_basin_boundaries');
                    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % potential
                elseif isfield(result, 'u') 
                    U = result.u(result.dom);
                    U = U - min(U) + 0.1 * (max(U) - min(U));
                    hax1 = fill([result.dom, result.dom(end), result.dom(1)], [U, 0, 0], [130/255 177/255 255/255]);
                    xlim(obj.domain);
                    xlabel(obj.namex);
                    ylabel('Potential')
                    if ~isfield(args, 'verticalbar') || args.verticalbar
                        obj.plot('vertical_basin_boundaries');
                    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % pdf                      
                elseif isfield(result, 'PDF') || isfield(result, 'pdfs')  
                    if isfield(args, 't')
                        t = args.t;
                    else
                        t = result.t(end);
                    end
                    hax1 = findobj(get(0, 'currentfigure'), 'tag', 'pdfs');
                    if ~isempty(hax1)
                        updatefig = true;
                    else
                        updatefig = false;
                    end
                    it = find(result.t >= t, 1);
                    if isempty(it)
                        it = length(result.t);
                    end
                    if (it == length(result.t)) && isfield(result, 'PDF')
                        if updatefig
                            set(hax1, 'YData', [result.PDF, 0 0]);
                        else
                            hax1 = fill([result.x result.x(end) result.x(1)], [result.PDF, 0 0], [130/255 177/255 255/255]);
                        end
                    else
                        if updatefig
                            set(hax1, 'YData', [result.pdfs(it, :), 0 0]);
                        else
                            hax1 = fill([result.x result.x(end) result.x(1)], [result.pdfs(it, :), 0 0], [130/255 177/255 255/255]);
                        end
                    end
                    xlim(obj.domain);
                    xlabel(obj.namex);
                    if isfield(result, 'PDF') && isempty(args)
                        ylabel('Stationary density')
                    else
                        ylabel('Probability density')
                    end
                    set(hax1, 'tag', 'pdfs');
                    if ~isfield(args, 'verticalbar') || args.verticalbar
                        obj.plot('vertical_basin_boundaries');
                    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % integration      
                elseif isfield(result, 'y')            
                    hax = plot(result.t, result.y, 'b-');
                    xlabel('Time');
                    ylabel(obj.namex);
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% % survival    
                elseif isfield(result, 'survival') 
                    [x, y] = meshgrid(result.x, result.t);
                    colormap( flipud([linspace(68, 243, 100)', linspace(138, 247, 100)', linspace(255, 255, 100)']./255));
                    [~, hax1] = contourf(x, y, result.survival, 100);
                    set(hax1, 'LineColor', 'none')
                    hold on
                    [~, hax1] = contour(x, y, result.survival, [0.5 0.5]);  % 
                    % [~, hax1] = contour(x, y, result.survival, [0.5 0.5], 'ShowText', 'on');
                    set(hax1, 'LineColor', 'k', 'LineWidth', 1.8)
                    [~, hax1] = contour(x, y, result.survival, [0.75 0.75]);
                    set(hax1, 'LineColor', [137 23 23]/256, 'LineWidth', 1.8)
                    [~, hax1] = contour(x, y, result.survival, [0.25 0.25]);
                    set(hax1, 'LineColor', [30 72 186]/256, 'LineWidth', 1.8)
                    % [~, hax1] = contour(x, y, result.survival, 1-[0.9:0.001:0.999]);
                    %  set(hax1, 'LineColor', 'g', 'LineWidth', 1)
                    xlim(obj.domain);
                    xlabel(strcat('Initial state',32, obj.namex));
                    % ylim auto
                    if isempty(obj.timeunit)
                        ylabel('Survival time')
                    else
                        ylabel(sprintf('Survival time (%s)', obj.timeunit))
                    end
                    if ~isfield(args, 'verticalbar') || args.verticalbar
                        obj.plot('vertical_basin_boundaries');
                    end
                elseif isfield(result, 'command')
                    obj.plot(result.command);
                end
            end
            if nargout > 0
                hax = hax1;
            end
        end
        
    end
end

% Euler integration
function [tout, yout] = euler(odefunct, tspan, y0, options)
    %Euler integration
    if (nargin < 4)
        %if there are no valid options use default stepsize
        delta = 0.01;
        options.StepSize = 0.3;
    else
        %use the option StepSize as step size
        if ~isfield(options, 'StepSize')
            options.StepSize = options.MaxStep;
        end
        if isempty(options.StepSize)
            options.StepSize = 0.1;
        end
        delta = options.StepSize;
    end
    if nargin < 4
        nonNegative = [];
    else
        nonNegative = odeget(options, 'NonNegative', [], 'fast'); 
    end
    anyNonNegative = ~isempty(nonNegative);
    haveOutputFcn = ~isempty(options.OutputFcn);
    outputFcn = options.OutputFcn;
    if haveOutputFcn
        feval(outputFcn, tspan, y0, 'init');
    end
    % Test that tspan is internally consistent.
    tspan = tspan(:);
    ntspan = length(tspan);
    if ntspan == 1
        t0 = 0;
        next = 1;
    else
        t0 = tspan(1);
        next = 2;
    end
    tfinal = tspan(ntspan);
    if t0 == tfinal
        error('euler:tpan', 'The last entry in tspan must be different from the first entry.');
    end
    tdir = sign(tfinal - t0);
    if any(tdir * (tspan(2:ntspan) - tspan(1:ntspan - 1)) <= 0)
        error('euler:tspan', 'The entries in tspan must strictly increase or decrease.');
    end
    t = t0;
    y = y0(:);
    neq = length(y);
    %adapt delta if there should be given more output
    step_tspan = median(diff(tspan));   
    delta = min(delta, step_tspan);
    
    % Set the output flag.
    outflag = ntspan > 2; % output only at tspan points

    % Allocate memory if we're generating output.
    delta = delta * tdir;

    if nargout > 0
        if outflag % output only at tspan points
            tout = tspan;
            outflag = delta ~= step_tspan;
            yout = nan(ntspan, neq);
        else
            tout = transpose(t0:delta:tfinal);
            if tout(end) < tfinal %if tfinal cannot divided in delta's
                tout(end + 1) = tfinal;
            end
            
            yout = nan(size(tout, 1), neq);
        end
        nout = 1;
        tout(nout) = t;
        yout(nout, :) = transpose(y);  
    end
    %
    %MAIN LOOP
    %evaluate the odefunction for the next time steps
    %fold=[];
    running = 1;
    while running
        f = feval(odefunct, t, y);  
        %  simple Euler
        ynew = y + f .* delta;
        if anyNonNegative
            ynew(nonNegative) = max(ynew(nonNegative), 0);
        end

        tnew = t + delta;
        if tnew >= tfinal
            running = 0;
        end
        
        if ~outflag % computed points, no refinement only the last value
            nout = nout + 1;
            if ~running
                if nout > length(tout)
                    nout = length(tout);
                end
                t1 = tout(nout);
                yout(nout, :) = transpose(y + (ynew - y) ./ (tnew - t) .* (t1 - t));
            else
                yout(nout, :) = transpose(ynew);
            end
            tnew = tout(nout);
            if haveOutputFcn
                if feval(outputFcn, tnew, transpose(ynew), '')
                    ndx = ~isnan(yout(:, 1));
                    yout = yout(ndx, :);
                    tout = tout(ndx);
                    running = false;
                end
            end
        elseif (tdir * (tnew - tspan(next)) >= 0) % at tspan, tspan assumed to be larger than delta
            nout = nout + 1;
            t1 = tout(nout);
            y1 = transpose((y + (ynew - y) ./ (tnew - t) .* (t1 - t)));
            yout(nout, :) = y1;
            next = next + 1;
            if haveOutputFcn
                if feval(outputFcn, t1, y, '')
                    ndx = ~isnan(yout(:, 1));
                    yout = yout(ndx, :);
                    tout = tout(ndx);
                    running = false;
                end
            end
        end
        
        y = ynew;
        t = tnew;
    end
end

function fun = integralx_fun(afun, ystart, xstart)   % h = integralx_fun(h1, obj.domain(1), obj.domain(2));   h1 = @(x)obj.D1(x)./obj.D2(x)
    if isa(afun, 'chebfun')
        fun = cumsum(afun);
    else
        if nargin < 2
            ystart = 0;
        end
        if nargin < 3
            xstart = 0;
        end
        h = @(x,t)afun(x);
        fun = @(timespan)myode45(h,timespan,ystart,xstart);
    end
end

function y = myode45(h, timespan, t0, x0)
    if length(timespan) == 1
        timespan = [t0 timespan];
        [~, y] = ode45(h, timespan, x0);
        y = y(end);
    else
        [~, y] = ode45(h, timespan(:), x0);
        if numel(y) < numel(timespan)   
            y(end + 1:numel(timespan)) = NaN;
        end
        y = reshape(y, size(timespan));
    end
end


function fun = diff_fun(afun, dt)
    if nargin < 3
        dt = 1E-7; %not too small for rounding errors
    end
    if isa(afun, 'chebfun')
        fun = diff(afun);  
    elseif isa(afun, 'function_handle')  
        fun = @(x)(afun(x+dt/2)-afun(x-dt/2))/dt;
    else
        error('langevin:diff', 'type of function not supported');
    end
end

function xout = simulate_vector(afun, x)
    xout = zeros(size(x));
    for i = 1:numel(x)
        xout(i) = afun(x(i));
    end
end


function rootsx = roots_fun(afun, dom)
    if isa(afun, 'chebfun')
        rootsx = roots(afun);
    else
        % find zeros without precision
        x = linspace(dom(1), dom(2), 100);
        res = afun(x);  % afun(x)为D1
        ndx = diff(sign(res)) ~= 0;   
        rootsx = x(ndx);   
        % adapt with precision
        for i = length(rootsx): -1:1
            if isnan(res(i))  
                rootsx(i) = [];
            else
                rootsx(i) = fzero(afun, rootsx(i));  
            end
        end
    end
end

% euler(odefun, tspan, y0, options)
function res = simulate_rng(randseed, odefun, tspan, y0, options)
    rng(randseed)  
    [res.t, res.y] = euler(odefun, tspan, y0, options);
end


% -------------------------------------------------------------------------
