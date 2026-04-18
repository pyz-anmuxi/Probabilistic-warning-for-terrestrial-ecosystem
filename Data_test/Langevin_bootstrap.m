function res = Langevin_bootstrap(data, reconstResult, bootmethod, nreplicates, boot_options)
    %Langevin_bootstrap - bootstrap  methods of LangevinReconst
    %This procedure implements different ways to assess the uncertainty of
    %LangevinReconstr. It includes block bootstrap and Monte-Carlo error propagation. 
    %
    %Usage: Langevin_bootstrap(data,reconstResult,bootmethod,nreplicates,
    %boot_options)
    %
    %data - original data set (can be empty)
    %reconstResult - structure with the results of LangevinReconstr (can
    %also be a langevin_eq with a real model, "Tau", "method", "domain",
    %"dt", "datasize" should then be in the boot_options)
    %bootmethod - string with kind of bootstrap method:
    %  'error-propagation' - use the reconstructed Langevin equation to
    %            generate replicate data sets for Monte Carlo error
    %            propagation - can assess bias (default method)
    %  'model-based' - use a model to generate the bootstrap data - can
    %            assess bias 
    %  'replicates' - the bootstrap data sets are provided in the boot_options 
    %  'block' - block bootstrap method (blocks of tau length) drawing with
    %            replication - simple percentile method for CL
    %nreplicates - number of replicate bootstrap data sets (default 100)
    %boot_options - struct with options or additional data. Fields:
    %   MaxStep - the stepsize for the Euler Maruyama solver (default 0.001)
    %   RemoveBias - remove bias from the results? (default true if possible) 
    %   model - langevin_eq object with the model for 'model-based' error
    %           propagation
    %   data - replicate data sets (columns with replicates) needed for
    %          'replicates' method
    %
    % by Egbert van Nes
    %
    % Latest version: https://git.wageningenur.nl/sparcs/langevinreconstruct
    %
    if nargin < 2
        error('LangevinRecostr:nodata', 'No enough arguments specified');
    end
    if nargin < 3
        bootmethod = 'error-propagation';
    end
    if nargin < 4
        nreplicates = 100;
    end
    if nargin < 5
        boot_options = [];
    end
    
    if isa(reconstResult, 'langevin_eq')
        modl = reconstResult;
        if ~isfield(boot_options, 'title')
            boot_options.title = '';
        end
        if ~isfield(boot_options, 'method')
            boot_options.method = 'Nadaraya-Watson';
        end
        if ~isfield(boot_options, 'bins')
            boot_options.bins = 50;
        end
        if ~isfield(boot_options, 'domain')
            boot_options.domain = modl.domain;
        end
        C = linspace(boot_options.domain(1), boot_options.domain(2), boot_options.bins);
        reconstResult = struct('D1', modl.D1(C), 'D2', modl.D2(C), 'D4', zeros(size(C)), 'ErrorD1', zeros(size(C)), 'ErrorD2', zeros(size(C)), 'N', zeros(size(C)), 'C', C, ...
            'options', struct('datasize', [1 max(boot_options.datasize)], 'domain', boot_options.domain, 'bins', boot_options.bins, 'Tau', boot_options.Tau, 'dt', boot_options.dt, 'method', boot_options.method, 'title', boot_options.title));
    end
    
    res = reconstResult;
    resbootstrap = struct('results', [], 'bootdata', [], 'bootmethod', bootmethod, 'nreplicates', nreplicates, 'options', boot_options);
    
    if iscell(bootmethod) 
        for i = 1:length(bootmethod)  
            resx = Langevin_bootstrap(data, reconstResult, bootmethod{i}, nreplicates, boot_options);
            resbootstrap(i) = resx.bootstrap;
        end
        res.bootstrap = resbootstrap;
        return;
    end
    
    L = reconstResult.options.domain(1);
    R = reconstResult.options.domain(2);
    bins = reconstResult.options.bins;
    Tau = reconstResult.options.Tau;
    dt = reconstResult.options.dt;
    method = reconstResult.options.method;
    
    if size(data, 2) < size(data, 1)
        data = data.';
    end
    
    if any(strcmpi(bootmethod, {'error-propagation', 'model-based', 'replicates'}))
        if ~isfield(boot_options, 'MaxStep')
            boot_options.MaxStep = 0.001;
        end
        if ~isfield(boot_options, 'RemoveBias')
            if strcmpi(bootmethod, 'error-propagation')
                boot_options.RemoveBias = true;
            else
                boot_options.RemoveBias = false;
            end
        end
       
        resbootstrap.options = boot_options;
        [res1, parentmod] = getresults(reconstResult, false);
        eqns = parentmod.equilibria;
        
        if strcmpi(bootmethod, 'model-based') || strcmpi(bootmethod, 'replicates')
            if ~(isfield(boot_options, 'model') || isfield(boot_options, 'data'))
                error('LangevinReconst_boot:nomodel', 'For model-based bootstrap, the option should have a field "model" containing a langevin_eq object\n or "data" with the data');
            end
            if isfield(boot_options, 'model')
                if ~isa(boot_options.model, 'langevin_eq')
                    error('LangevinReconst_boot:wrongmodel', 'The model should be a langevin_eq object');
                end
                parentmod = boot_options.model;
                eqns = parentmod.equilibria;
            else
                parentmod = [];
                bootstrapdatay = boot_options.data;
                reconstResult.options.data=[];
                resbootstrap.options.data=[];%otherwise the file becomes too big
               % boot_options=rmfield(boot_options,'data'); 
            end
        end
        
        if ~isempty(parentmod)
            Ndata = max(reconstResult.options.datasize);
            %draw initial conditions from the stationary pdf
            yinit = parentmod.draw_pdf(nreplicates * 2, 1);
            %create data
            fprintf('Creating %d bootstrap data sets of %d length\n', nreplicates * 2, Ndata);

            [bootstrapdata, fcn] = parentmod.simulate(linspace(0, Ndata * dt, Ndata), yinit, boot_options);
            bootstrapdatay = bootstrapdata.y;
           % resbootstrap.bootdata = bootstrapdatay;
            resbootstrap.bootdata=fcn; %we can reproduce the data but save much storage space
        else
       %     resbootstrap.bootdata = bootstrapdatay;
        end
        res1.haserror = true;
        % bootstrap(nreplicates) = res1;
        k = 0;
        tmpbootstrap = struct('reconstr', {}, 'results', {}, 'haserror', {});
        bootstrap = [];
        while length(bootstrap) < nreplicates && k < size(bootstrapdatay, 2) %two trials (this could be improved)
            n1 = nreplicates + floor(nreplicates / 10) - length(bootstrap); %a bit more to avoid short parfors
            tmpbootstrap(n1) = struct('reconstr', [], 'results', [], 'haserror', []);
            fprintf('n1 = %d\n', n1);
            remove_bias = boot_options.RemoveBias;
            for i = 1:n1
                if i + k <= size(bootstrapdatay, 2)
                    bootdata = bootstrapdatay(:, i + k);
                    try
                        res_boot = getresults(LangevinReconst(bootdata, L, R, bins, Tau, dt, method, sprintf('bootdata %d', i + k), false), remove_bias, eqns);
                        res_boot.haserror = false;
                        tmpbootstrap(i) = res_boot;
                    catch
                        %reject if there was an error calculating mean exit
                        %time
                        tmpbootstrap(i).haserror = true;
                    end
                else
                    tmpbootstrap(i).haserror = true;
                end
                disp(i);
            end
            bootstrap = [bootstrap, tmpbootstrap(~[tmpbootstrap(:).haserror])];
            fprintf('size bootstrap = %d\n', length(bootstrap));
            tmpbootstrap = struct('reconstr', {}, 'results', {}, 'haserror', {});
            k = k + n1;
        end
        if length(bootstrap) > nreplicates
            bootstrap = bootstrap(1:nreplicates);
        end
        resbootstrap.results = bootstrap;
        res.bootstrap = resbootstrap;
        %only save bootstraps that could find mean exit time
        if boot_options.RemoveBias
            if strcmpi(bootmethod, 'model-based') && isfield(boot_options, 'model')
                res = removebias(res, bootmethod, boot_options.model);
            else
                res = removebias(res);
            end
        end
    else
        
        if size(data, 1) == 1
            data = makesegments(data, max(Tau)).';
        end
        %get mean exit time with original data set 
        %  res1 = getresults(LangevinReconst(data, L, R, bins, Tau, dt, method));
        %res = reconstResult;
        %preallocate bootstrap data
        ndxes =  cell(nreplicates,1);
        %zeros(size(data, 2), nreplicates, 'int32');
        bootstrap(nreplicates) = struct('reconstr', [], 'results', [], 'haserror', []);
        parfor i = 1:nreplicates
            r1=rng;
            ndx= @()drawndx(r1,size(data,1),bootmethod);
            %ndx = drawndx(data, bootmethod);
            bootdata = data(ndx(),:);
            try
                %              res = getresults(LangevinReconst(bootdata, L, R, bins, Tau, dt, method, '', false));
                bootstrap(i) = getresults(LangevinReconst(bootdata, L, R, bins, Tau, dt, method, '', false), false);
            catch
                fprintf('Error in number %d, retry\nb', i);
                r1=rng;
                ndx= @()drawndx(r1,size(data,1),bootmethod);
                %ndx = drawndx(data, bootmethod);
                bootdata = data(ndx(),1);
                try
                    bootstrap(i) = getresults(LangevinReconst(bootdata, L, R, bins, Tau, dt, method, '', false), false);
                catch
                    bootstrap(i).haserror = true;
                end
            end
            ndxes{i} = ndx;
            disp(i);
        end
        resbootstrap.results = bootstrap(~[bootstrap(:).haserror]);
        res.bootstrap = resbootstrap;
        %res.bootstrap.bootdata = ndxes;
        %res.bootstrap.options.origdata = data;
    end
end
function [res1, mod] = getresults(res2, remove_bias, eqns)
     mod = langevin_eq(res2, 'weightedspline');  
    if nargin < 3
        eqns = mod.equilibria;
    end
    if ~remove_bias 
        pdf1 = mod.pdf;
        mod.mean_exit(eqns, pdf1);
    else
        mod.xtra.results = {};
    end
    res1 = mod.xtra; 
    res1.haserror = false;
end

function ndx = drawndx(randseed,n, bootmethod)
    rng(randseed)
    switch bootmethod
        case 'shuffle'
            %used for testing: only shuffle the data to see the effect
            %of long scale structure
            [~, ndx] = sort(rand(n, 1));
        otherwise
            %draw random segements with replacement
            ndx = randi(n, 1, n);
    end
end

function res = removebias(results, bootmethod, parentmodel)
    if ~isfield(results,'bootstrap')
        warning('Langinvin_boot:bias','Cannot remove bias if no bootstrap data are available');
        return;
    end
    if nargin < 2 && length(results.bootstrap) > 1
        bootmethod = 'error-propagation';
    end
    if nargin<3
        parentmodel=[];
    end
    if length(results.bootstrap) == 1
        bootmethod = results.bootstrap.bootmethod;
    end
    res = results;
    if isempty(parentmodel)
        parentmodel=langevin_eq(results,'weightedspline');
        eqns=parentmodel.equilibria;
        parent=results;
    elseif isfield(parentmodel.xtra,'reconstr')
        parent=parentmodel.xtra.reconstr;
        eqns=parentmodel.eqns;
    else
        parent=struct('D1',parentmodel.D1(results.C),'D2',parentmodel.D2(results.C),'D4',0);
        eqns=[];
    end
    % Fields of results.bootstrap(j)
    %        results: [1X100 struct]
    %       bootdata: [15001X200 double]
    %     bootmethod: 'model-based'
    %   nrepllicates: 100
    %        options: [1X1 struct]
    ndx = find(strcmp({results.bootstrap(:).bootmethod}, bootmethod));
    for j = ndx
        if isfield(res.bootstrap(j).options, 'bias_D1')
            warning('langevin:removebias:already', 'The bias was already removed');
        else
            boot_D1 = zeros(numel(results.bootstrap(j).results), numel(parent.D1));
            boot_D2 = zeros(numel(results.bootstrap(j).results), numel(parent.D1));
            boot_D4 = zeros(numel(results.bootstrap(j).results), numel(parent.D1));
            for i = 1:numel(results.bootstrap(j).results)
                boot_D1(i, :) = results.bootstrap(j).results(i).reconstr.D1;
                boot_D2(i, :) = results.bootstrap(j).results(i).reconstr.D2;
                boot_D4(i, :) = results.bootstrap(j).results(i).reconstr.D4;
            end
            res.bootstrap(j).options.bias.D1 = median(boot_D1, 1) - parent.D1; %negative is that the bootstrap median is lower than real D1
            res.bootstrap(j).options.bias.D2 = median(boot_D2, 1) - parent.D2;
            res.bootstrap(j).options.bias.D4 = median(boot_D4, 1) - parent.D4;
            res.bootstrap(j).options.parent=struct('D1',parent.D1,'D2',parent.D2,'D4',parent.D4);
            for i = 1:numel(results.bootstrap(j).results)
                res.bootstrap(j).results(i).reconstr.D1 = res.bootstrap(j).results(i).reconstr.D1 - res.bootstrap(j).options.bias.D1;
                res.bootstrap(j).results(i).reconstr.D2 = res.bootstrap(j).results(i).reconstr.D2 - res.bootstrap(j).options.bias.D2;
                res.bootstrap(j).results(i).reconstr.D4 = res.bootstrap(j).results(i).reconstr.D4 - res.bootstrap(j).options.bias.D4;
                if isfield(res.bootstrap(j).results(i),'results')
                    res.bootstrap(j).results(i).biased_results = res.bootstrap(j).results(i).results; %could not make this working in vector
                end
            end
            res1 = res.bootstrap(j).results;
            disp('Calculating unbiased exit times');
            parfor i = 1:numel(results.bootstrap(j).results)
                try
                    %              res = getresults(LangevinReconst(bootdata, L, R, bins, Tau, dt, method, '', false));
                    res2 = getresults_remove(res1(i).reconstr,eqns)
                    res1(i).results = res2.results;
                catch
                    res1(i).haserror = true;
                end
                disp(i);
            end
            for i = 1:numel(results.bootstrap(j).results)
                res.bootstrap(j).results(i).results = res1(i).results;
            end
        end
    end
end

function dat1 = makesegments(dat, len_segments)
    N = length(dat);
    dat1 = nan(len_segments + 1, N);
    dat1(1, :) = dat(1:N);
    for i = 1:len_segments
        dat1(i + 1, 1:N - i) = dat(1 + i:end);
    end
end

function [res1, mod] = getresults_remove(res2, eqns)
    mod = langevin_eq(res2, 'weightedspline');
    if nargin<2
        eqns=mod.equilibria;
    end
    pdf1 = mod.pdf;
    mod.mean_exit(eqns,pdf1);
    res1 = mod.xtra;
    res1.haserror = false;
end




