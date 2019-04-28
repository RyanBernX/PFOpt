function [A, b, x0, info] = gendatabd(varargin)
%EXPERIMENT.gendatabd  Generate data for a Blind Deconvolution experiment.
import util.*
ip = inputParser;
ip.addParameter('seed', 0);
ip.addParameter('channel', 'gray');
ip.addParameter('numcoeffs', 5003);
ip.addParameter('useoracle', true);
ip.addParameter('epsilon', 0);
ip.addParameter('normalize', true);
ip.addParameter('imagename', 'shapes-256');
ip.addParameter('filtername', 'motion');
ip.addParameter('fouriermeasurements', false);
ip.parse(varargin{:});

channel    = ip.Results.channel;
numcoeffs  = ip.Results.numcoeffs;
useoracle  = ip.Results.useoracle;
epsilon    = ip.Results.epsilon;
normalize  = ip.Results.normalize;
imagename  = ip.Results.imagename;
filtername = ip.Results.filtername;
fm         = ip.Results.fouriermeasurements;

rng(ip.Results.seed); % Set the RNG

switch lower(imagename)
   case 'shapes-256'
         image = imread('data/shapes.png');
   case 'goldballs-256'
         load data/goldballs-256.mat x0;
         image = abs(x0);
         clear x0;
   case 'goldballs-512'
         load data/goldballs-512.mat x0;
         image = abs(x0);
         clear x0;
   otherwise
         image = imread(['data/' imagename]);
      % error('[experiments.gendatabd] filtername ''%s'' invalid!',filtername);
end

switch lower(channel)
   case 'red'
      image = image(:,:,1);
   case 'green'
      image = image(:,:,2);
   case 'blue'
      image = image(:,:,3);
   case 'gray'
      image = mean(image,3);
end

image = double(image);

info.originalSignal = image;

if normalize
   image = image / normv(image);
else
   image = mat2gray(image);
end

switch lower(filtername)
   case 'motion'
      filter = kernel.motion(size(image),30*size(image,1)/256,45);
   case 'cmotion'
      filter = kernel.cmotion(size(image),11*size(image,1)/256);
   case 'dcmotion'
      filter = kernel.dcmotion(size(image),11*size(image,1)/256);
   otherwise
      error('[experiments.gendatabd] filtername ''%s'' invalid!',filtername);
end

info.originalKernel = filter;

if normalize
   filter = filter / normv(filter);
else
   filter = filter / sum(filter(:));
end

B2 = basis.matrix(filter>0);
x2 = B2.adjoint(filter);

H = basis.haar(true(size(image)));
c = H.adjoint(image);

xtmp = [c(:);conj(x2(:))];
Atmp = hop.bd(H,B2,isreal(xtmp),fm);
b0   = Atmp.forward(xtmp);

if ~useoracle
   if fm
      Atmp = hop.bd(H,B2,isreal(xtmp),false);
      btmp = Atmp.forward(xtmp);
      c    = H.adjoint(btmp);
   else
      c = H.adjoint(b0);
   end
end

[~,p] = sort(abs(c(:)),'descend');
support = false(size(image));
support(p(1:numcoeffs)) = true;
B1 = basis.haar(support);
x1 = B1.adjoint(image);

x0 = [x1(:);conj(x2(:))];

A = hop.bd(B1,B2,isreal(x0),fm);
b = A.forward(x0);

info.normalized = normalize;
info.epsilon = epsilon;
info.b0 = b;

end
