function stats = bd(varargin)
%EXPERIMENT.BD  Blind Deconvolution experiment.
%
% Options
%  solver: (saga)
%  solverOpts: ([]) options struct passed to the solver.
%  genOpts: ([]) options struct passed to the data generator.
%  plot: (false) plot the original and recovered solution.
import util.*
import kernel.*
p = inputParser;
p.addParameter('solver', 'saga');
p.addParameter('solverOpts', struct('iterations',1e3,'cutDim',2,...
                                    'optTol',1e-4,'feaTol',1e-4,...
                                    'LSmaxit',16,'LSeta',1,'skipDFP',false));
p.addParameter('genOpts', struct('normalize',true));
p.addParameter('plot', false);
p.parse(varargin{:});

genOpts = p.Results.genOpts;
solverOpts = p.Results.solverOpts;

[A, b, x0, info_gen] = experiments.gendatabd(genOpts);

global ndwt nidwt nfft2 nifft2;

switch p.Results.solver
   case 'saga'
      [x, r, info_sol] = saga_sd(A, b, solverOpts);
      signal = info_gen.originalSignal;
      kernel = info_gen.originalKernel;
      signalEst = A.b1.forward(x(A.i1))*normv(signal);
      kernelEst = A.b2.forward(conj(x(A.i2)))*normv(kernel);
      measurementsEst = A.forward(x)*normv(signal)*normv(kernel);
      x0 = [A.b1.adjoint(signal   );A.b2.adjoint(kernel   )];
      x  = [A.b1.adjoint(signalEst);A.b2.adjoint(kernelEst)];
      b  = A.forward(x0);
      r  = b-A.forward(x);

   case 'saga-feas'
      solverOpts.stopOnFeasible = true;
      solverOpts.skipDFP = true;
      [x, r, info_sol] = saga_sd(A, b, solverOpts);
      signal = info_gen.originalSignal;
      kernel = info_gen.originalKernel;
      signalEst = A.b1.forward(x(A.i1))*normv(signal);
      kernelEst = A.b2.forward(conj(x(A.i2)))*normv(kernel);
      measurementsEst = A.forward(x)*normv(signal)*normv(kernel);
      x0 = [A.b1.adjoint(signal   );A.b2.adjoint(kernel   )];
      x  = [A.b1.adjoint(signalEst);A.b2.adjoint(kernelEst)];
      b  = A.forward(x0);
      r  = b-A.forward(x);

   case 'aug-lag' % Ahmed, Recht, Romberg code
      [~,~,r,signal,kernel,measurements,signalEst,kernelEst,measurementsEst] = deblur(); close all;
      x0 = [A.b1.adjoint(signal   );A.b2.adjoint(fftshift(kernel   ))];
      x  = [A.b1.adjoint(signalEst);A.b2.adjoint(fftshift(kernelEst))];
      b  = A.forward(x0)*256;
      r  = b-A.forward(x)*256;
      nFFT = nfft2+nifft2;
      nDWT = ndwt+nidwt;
      info_sol.nfft = nFFT;
      info_sol.ndwt = nDWT;

end

% Compute errors
x01 =      x0(A.i1,:) ;
x02 = conj(x0(A.i2,:));
x1  =       x(A.i1,:) ;
x2  = conj( x(A.i2,:));

if ~info_gen.normalized
   x1 = x1*normv(x2,1)*sign(sum(real(x2(:))));
   x2 = x2/normv(x2,1)*sign(sum(real(x2(:))));
else
   phase = dotv(x01,x1)/normv(x01)/normv(x1);
   x1 = x1*phase;
   x2 = x2/phase;
end

xError = matrixerror(x1,conj(x2),x01,conj(x02),'fro');
xErrorRel = xError / normv(x01) / normv(x02); % = ||xx'-x0x0'||_F/||x0x0'||_F
x1Error = normv(x1-x01);
x1ErrorRel = x1Error / normv(x01); % = ||x1-x01||_F/||x01||_F
x2Error = normv(x2-x02);
x2ErrorRel = x2Error / normv(x02); % = ||x2-x02||_F/||x02||_F
rError = normv(r);
rErrorRel = rError / normv(b);

% Collect stats
stats.xError     = xError;
stats.xErrorRel  = xErrorRel;
stats.x1Error    = x1Error;
stats.x1ErrorRel = x1ErrorRel;
stats.x2Error    = x2Error;
stats.x2ErrorRel = x2ErrorRel;
stats.rError     = rError;
stats.rErrorRel  = rErrorRel;
stats.info_gen   = info_gen;
stats.info_sol   = info_sol;
stats.A          = A;
stats.x          = x;
stats.x0         = x0;
stats.x1         = x1;
stats.x2         = x2;
stats.x01        = x01;
stats.x02        = x02;
stats.r          = r;
stats.signalEst  = A.b1.forward(x1);
stats.kernelEst  = A.b2.forward(x2);
stats.bEst       = A.forward(x)*256;
stats.signal     = A.b1.forward(x01);
stats.kernel     = A.b2.forward(x02);
stats.b          = A.forward(x0)*256;

if p.Results.plot
   subplot(1,2,1); imshow(reshape(abs(A.b1.forward(x01)), A.m1, A.m2),[]);
   subplot(1,2,2); imshow(reshape(abs(A.b1.forward(x1 )), A.m1, A.m2),[]);
end

end
