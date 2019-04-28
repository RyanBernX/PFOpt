function [A, b, x0, info] = gendatapl(varargin)
%EXPERIMENT.gendatapl  Generate data for a PhaseLift experiment.
import util.*
ip = inputParser;
ip.addParameter('signal', 'gaussian');
ip.addParameter('mask', 'gaussian');
ip.addParameter('epsilon', 0);     % noise in data, factor of norm(b).
ip.addParameter('n', 128);
ip.addParameter('m', 1);
ip.addParameter('L', []);
ip.addParameter('normalize_signal', false); % set true to scale ||x0|| = 1.
ip.addParameter('normalize_mask', false); % set true to scale ||x0|| = 1.
ip.addParameter('regularize_mask', false);
ip.addParameter('seed', 0);
ip.addParameter('resizeImage', 1);
ip.addParameter('channel', 'gray');
ip.addParameter('imagefile', 'data/nebula.tif');
ip.parse(varargin{:});
n = ip.Results.n; % signal length
m = ip.Results.m; % signal width
L = ip.Results.L; % number of measurements
epsilon = ip.Results.epsilon;
mask = ip.Results.mask;
normalize_signal = ip.Results.normalize_signal;
normalize_mask = ip.Results.normalize_mask;
regularize_mask = ip.Results.regularize_mask;

rng(ip.Results.seed); % Set the RNG

% ---------------------------------------------------------------------
% Generate signal
% ---------------------------------------------------------------------
switch ip.Results.signal
   
   case '2Dreal'
      % This is the "Wusterland2" image from the Hubble Telescope:
      % http://hubblesite.org/newscenter/archive/releases/nebula/2015/12/image/a/
      % Use an octanary mask.
      x0 = imread(ip.Results.imagefile);
      if ip.Results.resizeImage < 1
         x0 = imresize(x0, ip.Results.resizeImage);
      end
      x0 = mat2gray(x0);
      switch lower(ip.Results.channel)
         case 'red'
            x0 = x0(:,:,1);
         case 'green'
            x0 = x0(:,:,2);
         case 'blue'
            x0 = x0(:,:,3);
         case 'gray'
            if size(x0, 3) == 3
               x0 = rgb2gray(x0);
            end
      end
      [n, m] = size(x0);
      mask = 'octanary';
      normalize_signal = false;
      
   case 'smooth'
      % Generate a smooth solution. This bit of code is taken from
      %    http://web.stanford.edu/~mahdisol/PRcode.html
      % It generates a random low-pass signal.
      m = 1;
      DFT_matrix = fftshift(fft(eye(n)))./sqrt(n);
      M = ceil(n/8);
      measurements_lowpass = (n/2-round(M/2-1)):(n/2+round(M/2));
      DFT_lowpass = DFT_matrix(measurements_lowpass,:);
      x0 = DFT_lowpass'*(randn(M,1) + 1i*randn(M,1));
      
   case {'gaussian','dual'}
      x0 = randn(n,m) + 1i*randn(n,m);

end

% Ensure signal is a vector (possibly normalized).
x0 = x0(:);
if normalize_signal
   x0 = x0 / norm(x0);
end

% ---------------------------------------------------------------------
% Generate mask.
% ---------------------------------------------------------------------
switch mask

   case 'binary'
      if isempty(L)
         L = 10;
      end
      M = cdp.make('binary', n, m, L, ...
         'forceregular',regularize_mask, 'normalize', normalize_mask);
      
   case 'gaussian'
      if isempty(L)
         L = 10;
      end
      M = cdp.make('cgaussian', n, m, L, ...
         'forceregular',regularize_mask, 'normalize', normalize_mask);

   case 'octanary'
      if isempty(L)
         L = 10;
      end
      M = cdp.make('octanary', n, m, L, ...
         'forceregular',regularize_mask, 'normalize', normalize_mask);
end

% ---------------------------------------------------------------------
% Create operator, measurements (RHS), etc.
% ---------------------------------------------------------------------
A = hop.pl(M,isreal(x0));  % phase-lift operator
b = A.forward(x0);         % measurements

if strcmp(ip.Results.signal,'dual')
   y = randn(size(b));
   [~,~,V,d] = A.gdobjective(y,1,struct('tol',eps,'p',5,'maxit',512));
   if d(1) < 0
      y = -y;
      [~,~,V,d] = A.gdobjective(y,1,struct('tol',eps,'p',5,'maxit',512));
   end
   if normalize_signal
      x0  = reshape(V(:,1),n,m);
   else
      x0  = reshape(V(:,1)/realsqrt(d(1)),n,m);
   end
   b   = A.forward(x0);
   eta = y/normv(y);
else
   eta = randn(size(b));
   eta = eta/normv(eta);
end

% Add noise.
b0 = b;
if epsilon > 0
   e       = epsilon;
   c       = [rdot(eta*(1+e),eta*(1-e)),...
                    -2*rdot(e*b0,e*eta),...
                       -rdot(e*b0,e*b0)];
   epsilon = max(real(roots(c)));
%    epsilon = e*normv(b0);
end
b = b0+epsilon*eta;

% Save stats to info
info.epsilon = epsilon;
info.b0 = b0;

end
