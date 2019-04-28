function [T, Tsumm] = noisypl(varargin)
%EXPERIMENT.NOISELESSPL  Noiseless PhaseLift experiments.
ip = inputParser;
ip.addParameter('verbosity', 0);
ip.addParameter('iterations', 5000);
ip.addParameter('test', false)
ip.addParameter('parallel', false);
ip.addParameter('signal', 'dual');
ip.addParameter('mask', 'octanary');
ip.addParameter('normalize_signal', true);
ip.addParameter('normalize_mask', true);
ip.parse(varargin{:});

signal           = ip.Results.signal;
mask             = ip.Results.mask;
normalize_signal = ip.Results.normalize_signal;
normalize_mask   = ip.Results.normalize_mask;

% Experimental parameters.
n = 128;
if ip.Results.test
   iter = 1;        % number of iterations per tuple
   Sseq = {'saga','wflow'}; % solver
   Lseq = {12,9,6};      % number of masks
   Eseq = {5e-1,5e-2,5e-3};   % epsilons
else
   iter = 100;                             % number of iterations per tuple
   Sseq = {'saga','saga-feas','wflow'};    % solver
   Lseq = {12,9,6};                        % number of masks
   Eseq = {5e-1,1e-1,5e-2,1e-2,5e-3,1e-3}; % epsilons
end
   
% Quiet solvers, max iterations.
solverOpts.verbosity  = ip.Results.verbosity;
solverOpts.iterations = ip.Results.iterations;
solverOpts.scale   = 1;
solverOpts.skipDFP = true;
solverOpts.optTol  = 2e-4;
solverOpts.feaTol  = 2e-4;
solverOpts.pfdopts = struct('gtol',2e-7,'ftol',0,'solver','minFunc');
% solverOpts.skipDFP = false;
% solverOpts.optTol  = 1e-4;
% solverOpts.feaTol  = 1e-4;
% solverOpts.pfdopts = struct('gtol',1e-8,'ftol',0,'solver','minFunc');
% solverOpts.dfpopts = struct('gtol',1e-4,'ftol',0);

% Reset RNG.
rng('default');

% Create an array of tuples: zip{i} = {solver, L, seed}
zip = util.cartproduct(Sseq,Lseq,Eseq,num2cell(1:iter));
ind = randperm(numel(zip));
zip = zip(ind);
nzip = length(zip);

% Create array to store stats.
data(nzip) = struct(...
   'solver',[],...
   'epsilon',[],...
   'L',[],...
   'seed',[],...
   'xError',[],...
   'xErrorRel',[],...
   'rError',[],...
   'rErrorRel',[],...
   'nfft',[],...
   'nAdjoint',[],...
   'nMeasure',[],...
   'time',[],...
   'status',[]);

% Disable the parfor if parallel is false.
if ip.Results.parallel
   p = parcluster;
   NumWorkers = p.NumWorkers;
else
   NumWorkers = 0;
end
parfor (i = 1:nzip, NumWorkers)
%for i=1:nzip
   [solver, L, epsilon, seed] = deal(zip{i}{:});

   genOpts = struct('n', n, 'L', L, 'epsilon', epsilon, 'seed', seed,...
                    'signal', signal, 'mask', mask, ...
                    'normalize_signal', normalize_signal,...
                    'normalize_mask', normalize_mask);
   
   %-------------------------------------------------------------------
   % Run experiment instance.
   %-------------------------------------------------------------------
   prefix = pwd;
   prefixCache = fullfile(prefix, 'cache');
   filename = fullfile(prefixCache, strrep(sprintf('solution_%s_%i_%0.3f_%i',zip{i}{:}),'.','_'));

   stats = run_get_experiment(filename, solver, solverOpts, genOpts, ...
      'solverOpts', 'genOpts');
   %-------------------------------------------------------------------

   data(i).solver = solver;
   data(i).epsilon = epsilon;
   data(i).L = L;
   data(i).seed = seed;
   data(i).xError = stats.xError;
   data(i).xErrorRel = stats.xErrorRel;
   data(i).rError = stats.rError;
   data(i).rErrorRel = stats.rErrorRel;
   data(i).nfft = stats.info_sol.nfft;
   data(i).nAdjoint = stats.info_sol.nAdjoint;
   data(i).nMeasure = stats.info_sol.nMeasure;
   data(i).time = stats.info_sol.time;
   data(i).status = stats.info_sol.status;
   
end % for

T = struct2table(data);
T.solver = categorical(T.solver);

succ = @(x)sum(x < 1e-2); % counts no. of successes

% Process data.
Tavg_nfft = cell(length(Sseq),1);
Tavg_xerr = cell(length(Sseq),1);
Tavg_succ = cell(length(Sseq),1);
for i = 1:length(Sseq)
    solver = Sseq{i};
    %if strcmp(solver,'saga'),  optmsg = 'Optimal solution found'; end
    %if strcmp(solver,'tfocs'), optmsg = 'Step size tolerance reached'; end
    solved = T.solver==solver;% & T.xErrorRel < 1e-2;% cellfun(@strcmp,T.status,repmat({optmsg},size(T.status)));
    Tavg_nfft{i}=varfun(@median,T(solved,:),'InputVariables',{'nfft'     },'GroupingVariables',{'solver','epsilon','L'});
    Tavg_xerr{i}=varfun(@median,T(solved,:),'InputVariables',{'xErrorRel'},'GroupingVariables',{'solver','epsilon','L'});
    Tavg_succ{i}=varfun( succ  ,T(solved,:),'InputVariables',{'xErrorRel'},'GroupingVariables',{'solver','epsilon','L'});
    % Sort by L, descending
    Tavg_nfft{i} = sortrows(sortrows(Tavg_nfft{i},'L','descend'),'epsilon');
    Tavg_xerr{i} = sortrows(sortrows(Tavg_xerr{i},'L','descend'),'epsilon');
    Tavg_succ{i} = sortrows(sortrows(Tavg_succ{i},'L','descend'),'epsilon');
end

% Print table
Lvals = Tavg_nfft{1}.L;
Evals = Tavg_nfft{1}.epsilon;
Tsumm = [];
for i = 1:length(Sseq)
   Tsumm = [Tsumm, ...
            Tavg_nfft{i}.median_nfft, ...
            Tavg_xerr{i}.median_xErrorRel, ...
            Tavg_succ{i}.Fun_xErrorRel  ]; %#ok<AGROW>
end
Tsumm = [Lvals Evals Tsumm];

fprintf(' %4s  %13s  %13s  %13s  %5s  %13s  %13s  %5s  %13s  %13s  %5s\n',...
   'L','eps','nFFT1','xErr1','succ1','nFFT2','xErr2','succ2','nFFT3','xErr3','succ3');
for i = 1:size(Tsumm,1)
   if ismember(Tsumm(i,1), [6 9 12])
      fprintf(' %4i  %13.4e  %13.4e  %13.4e  %5i  %13.4e  %13.4e  %5i  %13.4e  %13.4e  %5i\n',Tsumm(i,:)')
   end
end

T = sortrows(T,{'solver','epsilon','L'});

end % function 

function stats = run_get_experiment(filename, solver, solverOpts, genOpts, varargin)
   if exist([filename '.mat'],'file')
      saved_stats = load(filename);
      stats = saved_stats.stats;
      
   else
      stats = experiments.pl(...
         'solver', solver, ...
         'solverOpts', solverOpts, ...
         'genOpts', genOpts );
      
      % Remove r and x0 (waste of space).
      stats = rmfield(stats,{'r', 'x0'});

      save(filename, 'stats', varargin{:});
     
   end  
end
