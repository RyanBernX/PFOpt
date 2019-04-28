function [x, r, info] = saga_sd(A, b, varargin)

import util.*;

% save and reset global FFT counters
global nfft2 nifft2 ndwt nidwt;
if isempty(nfft2 ), nfft2  = 0; end
if isempty(nifft2), nifft2 = 0; end
if isempty(ndwt  ), ndwt   = 0; end
if isempty(nidwt ), nidwt  = 0; end
[dftPrevious,idftPrevious] = deal(nfft2,nifft2);
[dwtPrevious,idwtPrevious] = deal(ndwt ,nidwt );
[nfft2,nifft2,ndwt,nidwt]  = deal(0);

n = size(A, 2);
bNorm = normv(b);

% Parse input parameters.
p = inputParser;
p.addParameter('epsilon', 0);     % noise parameter
p.addParameter('iterations', n^2);
p.addParameter('maxStep', 1e+10);
p.addParameter('minStep', 1e-10);
p.addParameter('feaTol', 1e-5); % feasibility tol
p.addParameter('optTol', 1e-5); % optimality tol
p.addParameter('eigTol', 1e-2); % initial e-val tol for eigs
p.addParameter('eigTolFact', 0.1); % factor to reduct eigTol each itn
p.addParameter('eigIts', 300); % max number of iterations per eigs
p.addParameter('eigTolMin', eps); % min eigenvalue tol
p.addParameter('eigPlus', []);  % num of Lanczos vecs above cutDim
p.addParameter('eigInexact', 0);
p.addParameter('eigInexactDegree', 8);
p.addParameter('verbosity', 1);
p.addParameter('cutDim', 2);
p.addParameter('fid', 1);
p.addParameter('scale', []); % RHS of dual constraint: <b,y> = scale
p.addParameter('LSetamin', 0);  % Linesearch parameters
p.addParameter('LSetamax', 1);
p.addParameter('LSdelta', 1e-4);
p.addParameter('LSsigma', .9);
p.addParameter('LSrho', (sqrt(5)+1)/2);
% p.addParameter('LSeta', .85);
p.addParameter('LSeta', 1);
p.addParameter('LSmaxit', 16);   % 0 == no backtracks
p.addParameter('maxTime', 1800);% max number of seconds allowed for run
p.addParameter('StopOnFeasible', false);
p.addParameter('skipDFP', false);
p.addParameter('pfdOpts', struct()); % options for the PfD solver
p.addParameter('dfpOpts', struct()); % options for the DFP solver
p.parse(varargin{:});

epsilon = p.Results.epsilon;
fid = p.Results.fid;
verbosity = p.Results.verbosity;
maxItns = p.Results.iterations;
maxStep = p.Results.maxStep;
minStep = p.Results.minStep;
feaTol = p.Results.feaTol;
optTol = p.Results.optTol;
eigInexact = p.Results.eigInexact;
eigInexactDegree = p.Results.eigInexactDegree;
eigTol = p.Results.eigTol;
eigTolMin = p.Results.eigTolMin;
eigTolFact = p.Results.eigTolFact;
eigIts = p.Results.eigIts;
cutDim = p.Results.cutDim;
maxTime = p.Results.maxTime;
scale = p.Results.scale;
if isempty(scale)
   scale = bNorm - epsilon;
end
StopOnFeasible = p.Results.StopOnFeasible;
skipDFP = p.Results.skipDFP;
pfdOpts = p.Results.pfdOpts;
dfpOpts = p.Results.dfpOpts;

% Quiet subsolvers if saga_sd is quiet.
if verbosity >= 2
   pfdOpts.verbose = true;
   dfpOpts.verbose = true;
end

% Linesearch parameters.
LS.etamin = p.Results.LSetamin;
LS.etamax = p.Results.LSetamax;
LS.delta = p.Results.LSdelta;
LS.sigma = p.Results.LSsigma;
LS.rho = p.Results.LSrho;
LS.mu = p.Results.maxStep;
LS.eta = p.Results.LSeta;
LS.maxit = p.Results.LSmaxit;
LS.nls = 0;

% Exit conditions (constants).
stat = 0;
EXIT_OPTIMAL = 1;
EXIT_ITERATIONS = 2;
EXIT_INFEASIBLE_ITERATE = 3;
EXIT_TIME_OUT = 4;
EXIT_PRIMAL_INFEASIBLE = 5;
EXIT_OBJECTIVE_ERROR = 6;
EXIT_PRIMAL_FEASIBLE = 7;
EXIT_MSG = {
   'Unknown termination condition!'
   'Optimal solution found'
   'Too many iterations'
   'Infeasible dual iterate'
   'Time out'
   'Primal infeasible'
   'Objective error'
   'Feasible solution found (and requested)'
};

% Initialize counters.
nit = 0; % no. of iterations
nfe = 0; % no. of function evals

% Initialize eigenvalue options.
if isempty(p.Results.eigPlus)
   eigPlus = double(~A.isrealin);
else
   eigPlus = p.Results.eigPlus;
end
eigsOpts = struct('issym', true, ...
                  'isreal', A.isrealin, ...
                  'maxit', eigIts, ...
                  'p', min(n, max(2*cutDim+eigPlus, 20)), ...
                  'inexact', eigInexact, ...
                  'degree', eigInexactDegree);
time.objective = 0;
time.solver = 0;

% Start your engines!
tStart = tic;

%----------------------------------------------------------------------
% Initialize primal dual iterates, and associated gradient and step.
% 1. Construct a feasible point y, ie, <b, y> = scale.
% 2.
%----------------------------------------------------------------------
y = projection(b);
[duObj, g, v, lamDiff] = objective(y, A, [], eigsOpts, eigTol); nfe=nfe+1;
gNorm = normv(g, inf);

% Initialize linesearch quantities.
if gNorm < (1 / maxStep)
   step = maxStep;
else
   step = min( maxStep, max(minStep, 1/gNorm) );
end
LS.C = duObj;
LS.Q = 1;

%----------------------------------------------------------------------
% Log header.
%----------------------------------------------------------------------
logB = '%4i  %12.6e  %12.6e  %8.1e  %7.1e  %8.1e  %7.1e  %8.1e  %8.1e %3s %8.1e  %3i  %6i  %6i\n';
logH = '%4s  %12s  %12s  %8s  %7s  %7s  %8s  %8s  %8s %3s %8s  %3s  %6s  %6s\n';
if verbosity
   fprintf(fid,'\n');
   fprintf(fid,' %-20s: %6ix%6i %5s'   ,'matrix size'      ,n,n  ,'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'feasibility tol'  ,feaTol,'');
   fprintf(fid,' %-20s: %13i %5s'      ,'no. observations' ,numel(b),'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'optimality tol'   ,optTol   ,'');
   fprintf(fid,' %-20s: %13.2e %5s'    ,'max step'         ,maxStep,'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'scaling'          ,scale,'');
   fprintf(fid,' %-20s: %13.2e %5s'    ,'min step'         ,minStep,'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'norm b'           ,bNorm,'');
   fprintf(fid,' %-20s: %13i %5s'      ,'no. of eigs'      ,cutDim,'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'epsilon'          ,epsilon,'');
   fprintf(fid,' %-20s: %13i %5s'      ,'max its per eigs' ,eigIts,'');
   fprintf(fid,' %-20s: %13.2e %5s\n'  ,'eigTol decrease fact',eigTolFact,'');
   fprintf(fid,'\n');
end

%----------------------------------------------------------------------
% MAIN LOOP.
%----------------------------------------------------------------------

% Primal recovery.
[x, r, pfdstats] = primal_recovery(A, b, epsilon, y, v, g);

% Improve dual estimate.
spaced = '';
if skipDFP
   dfpstats.fOut = 0;
   dfpstats.iterations = 0;
else
   [yT, duObjT, gT, vT, lamDiffT, dfpstats] = dual_recovery(A, b, x, y, v);
   if duObjT < LS.C
      LS.C    = (duObjT-duObj)/LS.Q+LS.C;
      y       = yT;
      duObj   = duObjT;
      g       = gT;
      v       = vT;
      lamDiff = lamDiffT;
      spaced  = 'Y';
   end
end

while true

   % Quantities needed for log and to test exit conditions.
   duFeas = log(min(1, rdot(y,b) - epsilon*normv(y)));
   prFeas = max(0,normv(r)-epsilon);
   prObj = normv(x)^2;
   mGap = abs(1-prObj*duObj/scale);

   testPrFeas = prFeas <= feaTol*(1+normv(b));
   testOptiml = abs(mGap) <= optTol;

   if ~stat  &&  duObj < 0
      stat = EXIT_PRIMAL_INFEASIBLE;
   end

   if ~stat  &&  duFeas > feaTol
      stat = EXIT_INFEASIBLE_ITERATE;
   end

   if ~stat  && prFeas  &&  StopOnFeasible
      stat = EXIT_PRIMAL_FEASIBLE;
   end

   if ~stat  &&  testPrFeas && testOptiml
      stat = EXIT_OPTIMAL;
   end

   if ~stat  &&  nit >= maxItns
      stat = EXIT_ITERATIONS;
   end

   if ~stat  &&  toc(tStart) >= maxTime
      stat = EXIT_TIME_OUT;
   end

   %------------------------------------------------------------------
   % Print log and act on exit conditions.
   %------------------------------------------------------------------
   if verbosity
      if mod(nit, 1) == 0
         fprintf(fid,'\n');
         fprintf(fid,logH,...
            'itn','1/pr Obj','du Obj','mGap','l1-l2','step',...
            'prFeas','PFDpgNrm','dfp','spc','eigTol','nls',...
            'PFDits','DFPits');
      end
      fprintf(fid, logB, ...
         nit, 1/prObj, duObj/scale, mGap, lamDiff, step, prFeas,...
         pfdstats.pgNorms(end), realsqrt(2*dfpstats.fOut), ...
         spaced, eigTol, LS.nls, pfdstats.iterations,...
         dfpstats.iterations);
   end

   if stat
      break
   end

   %==================================================================
   % New iteration starts here.
   %==================================================================
   nit = nit + 1;

   %------------------------------------------------------------------
   % Linesearch.
   %------------------------------------------------------------------
   eigTol = max(eigTolMin, min(eigTol, eigTolFact*min([prFeas, mGap, eigTol])));
   yk = y; gk = g;
   try
      [duObj, y, g, v, lamDiff, LS] = linesearch(@objective, @projection,...
         yk, gk, step, A, v, eigsOpts, eigTol, LS);
   catch
      stat = EXIT_OBJECTIVE_ERROR;
   end
   nfe = nfe + LS.nls;

   % Primal recovery.
   [x, r, pfdstats] = primal_recovery(A, b, epsilon, y, v, g);

   % Spacer step.
   spaced = '';
   if ~skipDFP
      [yT, duObjT, gT, vT, lamDiffT, dfpstats] = dual_recovery(A, b, x, y, v);
      if duObjT < LS.C
         LS.C    = (duObjT-duObj)/LS.Q+LS.C;
         y       = yT;
         duObj   = duObjT;
         g       = gT;
         v       = vT;
         lamDiff = lamDiffT;
         spaced  = 'Y';
      end
   end

   % Update Barzilai-Borwein step length.
   dy = y - yk;
   dg = g - gk;
   step = min(maxStep, max(minStep, rdot(dy,dy)/rdot(dy,dg)));

end % while true
time.solver = toc(tStart);

%----------------------------------------------------------------------
% Log footer.
%----------------------------------------------------------------------
nFFT = nfft2+nifft2;
nDWT = ndwt +nidwt ;

if verbosity
fprintf(fid,'\n EXIT -- %s\n\n',EXIT_MSG{stat+1});
fprintf(fid,' %-22s: %10i %5s'  ,'no. of measurements', A.nForward,'');
fprintf(fid,' %-22s: %10.2f %5s\n','time solve (sec)', time.solver,'');
fprintf(fid,' %-22s: %10i %5s'  ,'no. of adjoints', A.nAdjoint,'');
fprintf(fid,' %-22s: %10.2f %5s\n','time objective (sec)', time.objective,'');
fprintf(fid,' %-22s: %10i %5s'  ,'no. of FFTs', nFFT,'');
fprintf(fid,' %-22s: %10i %5s\n','no. of PFDObj. calls', A.nPFDObjective,'');
fprintf(fid,' %-22s: %10i %5s'  ,'no. of DWTs', nDWT,'');
fprintf(fid,' %-22s: %10i %5s\n','no. of DFPObj. calls', A.nDFPObjective,'');
fprintf(fid,'\n');
end

% Gather exit info.
info.prFeas = prFeas;
info.duObj = duObj;
info.prObj = prObj;
info.mGap = mGap;
info.epsilon = epsilon;
info.scale = scale;
info.nit = nit;
info.nfe = nfe;
info.nfft = nFFT;
info.ndwt = nDWT;
info.nAdjoint = A.nAdjoint;
info.nMeasure = A.nForward;
info.time = time.solver;
info.stat = stat;
info.status = EXIT_MSG{stat+1};
info.y = y;

% restore global FFT and DWT counters
[nfft2,nifft2,ndwt,nidwt] = deal(dftPrevious,idftPrevious,dwtPrevious,idwtPrevious);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NESTED FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function [fd, g, v, lamDiff] = objective(y, A, v, eigsOpts, eTol)
      %OBJECTIVE  The dual objective value lambda_1(A'*y).
      %
      % JUL 6 2015: EIGS routine can fail if the parameter `p` is too
      % small. The catch statment below will detect a failure and
      % increase `p` 3 times in a row. Otherwise, will throw the error.

      tObjStart    = tic;
      eigsOpts.tol = eTol;
      eTries = 0; % No. of times eigenvalue solver called.
      while true
         try
            [fd,g,V,d] = A.gdobjective(y,cutDim,eigsOpts,v);
            break
         catch ME
            if eTries <= 3 && strcmp(ME.identifier,'MATLAB:eigs:ARPACKroutineErrorMinus14')
               eigsOpts.p = 2*eigsOpts.p;
               eTries = eTries + 1;
            else
               causeException = MException('MATLAB:saga:objective','eigs failed');
               ME = addCause(ME,causeException);
               rethrow(ME)
            end
         end
      end
      v = V(:,1);
      if numel(d) > 1
         lamDiff = (d(1)-d(2))/d(1);
      else
         lamDiff = NaN;
      end
      time.objective = time.objective + toc(tObjStart);
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function z = projection(z)
      % Projection operator onto the constraint
      % <b,y> - epsilon||y|| >= scale.
      z = project(z, b, epsilon, scale);
%     z = z + ((scale-rdot(b,z))/bNorm^2)*b; % projection onto <b,y>=1.
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function [x, r, pfdstats] = primal_recovery(A, b, epsilon, y, v, g)
      %PRIMAL_REOCOVERY  Recovery primal estimate from gradient.
      [x, pfdstats] = pfd(A, b, epsilon, y, v, g, pfdOpts);
      r = b - A.forward(x);
   end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   function [y, duObj, g, v, lamDiff, dfpstats] = dual_recovery(A, b, x, y, v)
      %DUAL_RECOVERY  Recover an dual estimate from primal.
      normr = util.normv(r);
      if epsilon > 0 && normr > 0
         y0 = r*(scale/(util.rdot(r,b)-epsilon*normr));
      else
         y0 = y;
      end
      [y, dfpstats] = dfp(A, b, epsilon, scale, x, y0, dfpOpts);
      [duObj, g, v, lamDiff] = objective(y, A, v, eigsOpts, eps);
   end

end % function saga_sd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCAL FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [duObj, y, g, v, lamDiff, LS] = linesearch(objective, projection, yk, gk, step, A, v, eigsOpts, eigTol, LS)
   % Full nonmonote linesearch.
   import util.*;
   y = projection( yk - step * gk );
   [duObj, g, v, lamDiff] = objective(y, A, v, eigsOpts, eigTol);
   dy = y - yk;
   done = false;
   if step <= LS.mu && duObj <= LS.C + LS.delta * rdot(gk,dy)
      % Armijo satisfied
      if rdot(g,dy) >= LS.sigma * rdot(gk,dy)
         % Wolfe satisfied
         done = true;
      else
         increaseStep = true;
      end
   else
      increaseStep = false;
   end
   nls = 1;
   while ~done && nls <= LS.maxit
      if increaseStep
         step = step * LS.rho;
         yT = projection( yk - step * gk );
         [duObjT, gT, vT, lamDiffT] = objective(yT, A, v, eigsOpts, eigTol);
         dyT = yT - yk;
         armijoSatisfied = (step <= LS.mu) && (duObjT <= LS.C + LS.delta * rdot(gk,dyT));
         done = ~armijoSatisfied || (rdot(gT,dyT) >= LS.sigma * rdot(gk,dyT));
         if armijoSatisfied || nls == LS.maxit
            y = yT; duObj = duObjT; g = gT; v = vT; lamDiff = lamDiffT;
         end
      else
         step = step / LS.rho;
         y = projection( yk - step * gk );
         [duObj, g, v, lamDiff] = objective(y, A, v, eigsOpts, eigTol);
         dy = y - yk;
         done = (duObj <= LS.C + LS.delta * rdot(gk,dy));
      end
      nls = nls + 1;
   end
   % Update linesearch quantities.
   LS.nls = nls;
   LS.Q = LS.eta * LS.Q + 1;
   LS.C = ( duObj + LS.C * (LS.Q - 1) ) / LS.Q;
end

