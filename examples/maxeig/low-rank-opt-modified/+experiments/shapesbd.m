function shapesbd(varargin)
%EXPERIMENT.SHAPESBD  Blind deconvolution on the shapes problem.
ip = inputParser;
ip.addParameter('verbosity', 0);
ip.addParameter('iterations', 1000);
ip.addParameter('test', false);
ip.addParameter('parallel', false);
ip.addParameter('plot', false);
ip.parse(varargin{:});

% Experimental parameters.
if ip.Results.test
   Sseq = {'saga','saga-feas'}; % solver
else
   Sseq = {'saga','saga-feas','aug-lag'}; % solver
end
   
% Quiet solvers, max iterations.
solverOpts = struct(...
   'verbosity', ip.Results.verbosity,...
   'iterations', ip.Results.iterations,...
   'optTol', 5e-2, 'feaTol', 1e-5,...
   'cutDim', 1, 'eigTol', eps, 'eigPlus', 0, 'LSmaxit', 0,...
   'pfdOpts', struct('gTol',1e-9), 'dfpOpts', struct('gTol',1e-7));

% Reset RNG.
rng('default');

% Create an array of tuples: zip{i} = {solver, L, seed}
zip = util.cartproduct(Sseq);
ind = randperm(numel(zip));
zip = zip(ind);
nzip = length(zip);

% Create array to store stats.
data(nzip) = struct(...
   'solver',[],...
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
%parfor (i = 1:nzip, NumWorkers)
for i=1:nzip
   [solver] = deal(zip{i}{:});

   genOpts = struct();
                 
   %-------------------------------------------------------------------
   % Run experiment instance.
   %-------------------------------------------------------------------
   prefix = pwd;
   prefixCache = fullfile(prefix, 'cache');
   filename = fullfile(prefixCache, sprintf('shapesbd_%s',zip{i}{:}));

   stats = run_get_experiment(filename, solver, solverOpts, genOpts, ...
      'solverOpts', 'genOpts');
   %-------------------------------------------------------------------

   data(i).solver = solver;
   data(i).xError = stats.xError;
   data(i).xErrorRel = stats.xErrorRel;
   data(i).x1Error = stats.x1Error;
   data(i).x1ErrorRel = stats.x1ErrorRel;
   data(i).x2Error = stats.x2Error;
   data(i).x2ErrorRel = stats.x2ErrorRel;
   data(i).rError = stats.rError;
   data(i).rErrorRel = stats.rErrorRel;
   data(i).nfft = stats.info_sol.nfft;
   data(i).ndwt = stats.info_sol.ndwt;
   data(i).stats = stats;

   indices = fftshift(reshape(1:2^16,256,256));
   stats.kernel    = stats.kernel(indices);
   stats.kernelEst = stats.kernelEst(indices);

   imwrite(uint8(round(max(0,min(255,stats.signal   ))    )),sprintf('cache/%s_signal.png'   ,solver));
   imwrite(uint8(round(255*stats.kernel/max(stats.kernel(:)))),sprintf('cache/%s_kernel.png'   ,solver));
   imwrite(uint8(round(max(0,min(255,stats.b        ))    )),sprintf('cache/%s_measur.png'   ,solver));
   imwrite(uint8(round(max(0,min(255,stats.signalEst))    )),sprintf('cache/%s_signalEst.png',solver));
   imwrite(uint8(round(255*stats.kernelEst/max(stats.kernel(:)))),sprintf('cache/%s_kernelEst.png',solver));
   imwrite(uint8(round(max(0,min(255,stats.bEst     ))    )),sprintf('cache/%s_measurEst.png',solver));
end % for

% Solver  nDFT  nDWT  x1Err  x2Err  resid
logB =  '%10s  %5i  %5i  %10.2e  %10.2e  %10.2e\n';
fprintf('%10s  %5s  %5s  %10s  %10s  %10s\n',...
        'solver','nDFT','nDWT','x1Err','x2Err','resid');
for i=1:length(Sseq)
    fprintf(logB, data(i).solver, data(i).nfft, data(i).ndwt, ...
                  data(i).x1ErrorRel, data(i).x2ErrorRel, ...
                  data(i).rErrorRel);
end

end % function 

function stats = run_get_experiment(filename, solver, solverOpts, genOpts, varargin)
   if exist([filename '.mat'],'file')
      saved_stats = load(filename);
      stats = saved_stats.stats;
      
   else
      stats = experiments.bd(...
         'solver', solver, ...
         'solverOpts', solverOpts, ...
         'genOpts', genOpts );
      
      % Remove r and x0 (waste of space).
      stats = rmfield(stats,{'r', 'x0'});
      
      save(filename, 'stats', varargin{:});
     
   end  
end
