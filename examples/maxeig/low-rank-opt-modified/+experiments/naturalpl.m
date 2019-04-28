function [T, Tsumm] = naturalpl(varargin)
%EXPERIMENT.NATURALPL  Natural image PhaseLift experiment.
ip = inputParser;
ip.addParameter('verbosity', 0);
ip.addParameter('iterations', 1000);
ip.addParameter('test', false)
ip.addParameter('parallel', false);
ip.addParameter('maxTime',120*60);
ip.addParameter('resizeImage',1);
ip.parse(varargin{:});
resizeImage = ip.Results.resizeImage;

% Experimental parameters.
if ip.Results.test
   Sseq = {'saga','saga-feas','wflow'}; % solver
   Lseq = {10, 15};                     % number of masks
   channel = 'gray';
   resizeImage = .1;                    % smaller test image
else
   Sseq = {'saga','saga-feas','wflow'}; % solver
   Lseq = {10, 15};                     % number of masks
   channel = 'gray';
 end
   
% Reset RNG.
rng('default');

% Create an array of tuples: zip{i} = {solver, L, seed}
zip = util.cartproduct(Sseq,Lseq);
ind = randperm(numel(zip));
zip = zip(ind);
nzip = length(zip);

% Create array to store stats.
data(nzip) = struct(...
   'solver',[],...
   'L',[],...
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
   p = parpool('local',3);
   NumWorkers = p.NumWorkers;
   warning off
   pctRunOnAll maxNumCompThreads(8)
   warning off
else
   NumWorkers = 0;
end
%parfor (i = 1:nzip, NumWorkers)
for i=1:nzip
   [solver, L] = deal(zip{i}{:});

   genOpts = struct('signal','wusterland','L', L, 'seed', L, ...
                    'channel',channel,'resizeImage',resizeImage);

   %-------------------------------------------------------------------
   % Solver-specific options.
   %-------------------------------------------------------------------
   solverOpts = struct('verbosity',ip.Results.verbosity,...
                       'iterations',ip.Results.iterations,...
                       'maxTime',ip.Results.maxTime);
   if strcmpi(solver,'wflow')
      solverOpts.gTol = 1e-7;
   end

   %-------------------------------------------------------------------
   % Run experiment instance.
   %-------------------------------------------------------------------
   prefix = pwd;
   prefixCache = fullfile(prefix, 'cache');
   filename = fullfile(prefixCache, sprintf('wusterland_%s_%i',zip{i}{:}));

   stats = run_get_experiment(filename, solver, solverOpts, genOpts, ...
      'solverOpts', 'genOpts');
   %-------------------------------------------------------------------

   data(i).solver = solver;
   data(i).L = L;
   data(i).rGap = stats.info_sol.mGap / max(1,norm(stats.info_gen.b0(:)));
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

% Process data.
Tavg_nfft = cell(length(Sseq),1);
Tavg_xerr = cell(length(Sseq),1);
Tavg_rgap = cell(length(Sseq),1);
for i = 1:length(Sseq)
    solver = Sseq{i};
    solved = T.solver==solver;
    Tavg_nfft{i}=varfun(@median,T(solved,:),'InputVariables',{'nfft'     },'GroupingVariables',{'solver','L'});
    Tavg_xerr{i}=varfun(@median,T(solved,:),'InputVariables',{'xErrorRel'},'GroupingVariables',{'solver','L'});
    Tavg_rgap{i}=varfun(@median,T(solved,:),'InputVariables',{'rGap'     },'GroupingVariables',{'solver','L'});
    % Sort by L, descending
    Tavg_nfft{i} = sortrows(Tavg_nfft{i},'L','descend');
    Tavg_xerr{i} = sortrows(Tavg_xerr{i},'L','descend');
    Tavg_rgap{i} = sortrows(Tavg_rgap{i},'L','descend');
end

% Print table
Lvals = Tavg_nfft{1}.L;
Tsumm = [];
for i = 1:length(Sseq)
   Tsumm = [Tsumm, ...
            Tavg_nfft{i}.median_nfft, ...
            Tavg_xerr{i}.median_xErrorRel, ...
            Tavg_rgap{i}.median_rGap  ]; %#ok<AGROW>
end
Tsumm = [Lvals Tsumm];

fprintf(' %4s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s\n',...
   'L','nFFT1','xErr1','rGap1','nFFT2','xErr2','rGap2','nFFT3','xErr3','rGap3');
for i=1:size(Tsumm,1)
   fprintf(' %4i  %14.5e  %14.5e  %14.5e  %14.5e  %14.5e  %14.5e  %14.5e  %14.5e  %14.5e\n',Tsumm(i,:)')
end

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
