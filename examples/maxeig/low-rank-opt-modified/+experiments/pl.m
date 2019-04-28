function stats = pl(varargin)
%EXPERIMENT.PL  PhaseLift experiment.
%
% Options
%  solver: (saga)
%  solverOpts: ([]) options struct passed to the solver.
%  genOpts: ([]) options struct passed to the data generator.
%  plot: (false) plot the original and recovered solution.
import util.*
p = inputParser;
p.addParameter('solver', 'saga');
p.addParameter('solverOpts', struct());
p.addParameter('genOpts', struct());
p.addParameter('plot', false);
p.parse(varargin{:});

genOpts = p.Results.genOpts;
solverOpts = p.Results.solverOpts;

[A, b, x0, info_gen] = experiments.gendatapl(genOpts);

if ~isfield(solverOpts,'epsilon') || isempty(solverOpts.epsilon)
   solverOpts.epsilon = info_gen.epsilon;
end

switch p.Results.solver
   case 'saga'
      [x, r, info_sol] = saga_sd(A, b, solverOpts);

   case 'saga-feas'
      solverOpts.stopOnFeasible = true;
      solverOpts.cutDim = 1;
      solverOpts.eigPlus = 1;
      [x, r, info_sol] = saga_sd(A, b, solverOpts);

   case 'tfocs'
      [x, r, info_sol] = tfocs_pl(A, b, solverOpts);
      info_sol.mGap = nan; % not reported by tfocs
      
   case 'wflow'
      % The `wflow` routine does its own measuring. Otherwise, we have
      % to redefine its internal operators, which isn't advisable!
      X0 = reshape(x0, A.n1, A.n2);
      if solverOpts.epsilon > 0
         solverOpts.Y = b*A.n;
      end
      [X, info_sol] = wflow(conj(A.masks), X0, solverOpts);
      r = b - A.forward(X);
      x = X(:);
      info_sol.mGap = nan; % not relevant for wflow
end

% Compute errors
xError = hermitianerror(x(:),x0(:),'fro');
xErrorRel = xError / normv(x0)^2; % = ||xx'-x0x0'||_F/||x0x0'||_F
rError = normv(r);
rErrorRel = rError / normv(b);

% Collect stats
stats.xError = xError;
stats.xErrorRel = xErrorRel;
stats.rError = rError;
stats.rErrorRel = rErrorRel;
stats.info_gen = info_gen;
stats.info_sol = info_sol;
stats.x = x;
stats.x0 = x0;
stats.r = r;

if p.Results.plot
   subplot(1,2,1); imshow(reshape(abs(x0), A.n1, A.n2));
   subplot(1,2,2); imshow(reshape(abs(x ), A.n1, A.n2));
end

end
