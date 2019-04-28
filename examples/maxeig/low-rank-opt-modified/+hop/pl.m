classdef pl < hop.hop

   properties
      masks
      m1
      m2
      m3
      m
      n1
      n2
      n
      eigWork
   end % properties

   methods

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function self = pl(masks,isrealin)
         if nargin < 2
            isrealin = false;
         end
         self@hop.hop(isrealin);
         self.masks = masks;
         self.m1    = size(self.masks,1);
         self.m2    = size(self.masks,2);
         self.m3    = size(self.masks,3);
         self.m     = self.m1*self.m2*self.m3;
         self.n1    = self.m1;
         self.n2    = self.m2;
         self.n     = self.n1*self.n2;

         self.eigWork = struct;
         self.eigWork.force_exact = 5;
      end % pl (constructor)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function out = forward(self,x1,x2,scatter)
         %FORWARD  Forward PhaseLift operator.
         if (nargin < 4) || isempty(scatter)
            scatter = false;
         end
         x1 = reshape(x1,[self.n1,self.n2,1,numel(x1)/self.n]);
         X1 = fft2(bsxfun(@times,self.masks,x1))/realsqrt(self.n);
         util.addToNFFT2(numel(X1)/self.n);
         if (nargin >= 3) && ~isempty(x2)
            x2 = reshape(x2,[self.n1,self.n2,1,numel(x2)/self.n]);
            if scatter
               x2 = permute(x2,[1,2,3,5,4]);
            end
            X2 = fft2(bsxfun(@times,self.masks,x2))/realsqrt(self.n);
            util.addToNFFT2(numel(X2)/self.n);
         else
            X2 = X1;
         end
         out = bsxfun(@times,X1,conj(X2));
         if ~scatter
            out = sum(out,4);
         end
         self.nForward = self.nForward + 1;
      end % forward

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function out = adjoint(self,y,x)
         x   = reshape(x,[self.n1,self.n2,      1,numel(x)/self.n,              1]);
         y   = reshape(y,[self.m1,self.m2,self.m3,              1,numel(y)/self.m]);
         out = sum(bsxfun(@times,conj(self.masks),ifft2(bsxfun(@times,y,fft2(bsxfun(@times,self.masks,x))))),3);
         out = reshape(out,[self.n1,self.n2,numel(x)/self.n,numel(y)/self.m]);
         util.addToNFFT2 (self.m3*numel( x )/self.n);
         util.addToNIFFT2(self.m3*numel(out)/self.n);
         self.nAdjoint = self.nAdjoint + 1;
      end % adjoint

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [f,g,V,d] = gdobjective(self,y,k,opts,v)
         if (nargin <= 4) || (numel(v) ~= self.n) || (isreal(v) ~= self.isrealin)
            v = [];
         end
         if nargin <= 3
            opts = [];
         end
         if nargin <= 2
            k = [];
         end
         if ~isempty(v)
            eigsopts.v = v;
         elseif ~isempty(opts) && isfield(opts,'v') && (numel(opts.v) == self.n) && (isreal(opts.v) == self.isrealin)
            eigsopts.v = opts.v;
         end
         if isempty(k)
            if ~isempty(opts) && isfield(opts,'cutDim') && ~isempty(opts.cutDim) && (opts.cutDim >= 1)
               k = opts.cutDim;
            else
               k = 1;
            end
         end
         if ~isempty(opts) && isfield(opts,'eigTol') && ~isempty(opts.eigTol)
            eigsopts.tol = opts.eigTol;
         end
         if ~isempty(opts) && isfield(opts,'eigIts') && ~isempty(opts.eigIts)
            eigsopts.maxit = opts.eigIts;
         end
         if ~isempty(opts) && isfield(opts,'eigPlus') && ~isempty(opts.eigPlus)
            eigsopts.p = 2*k+opts.eigPlus;
         end
         if ~isempty(opts) && isfield(opts,'tol') && ~isempty(opts.tol)
            eigsopts.tol = opts.tol;
         end
         if ~isempty(opts) && isfield(opts,'maxit') && ~isempty(opts.maxit)
            eigsopts.maxit = opts.maxit;
         end
         if ~isempty(opts) && isfield(opts,'p') && ~isempty(opts.p)
            eigsopts.p = opts.p;
         end
         if ~isempty(opts) && isfield(opts,'inexact') && ~isempty(opts.inexact)
             inexact = opts.inexact;
         else
             inexact = 0;
         end
         if ~isempty(opts) && isfield(opts,'degree') && ~isempty(opts.degree)
             self.eigWork.degree = opts.degree;
         else
             self.eigWork.degree = 8;
         end
         if self.isrealin
            Afun  = @(x) reshape(util.rvec(self.adjoint(y,x)), size(x));
            sigma = 'LA';
         else
            Afun  = @(x) reshape(util.vec(self.adjoint(y,x)), size(x));
            sigma = 'LR';
         end

         eigsopts.issym  = true;
         eigsopts.isreal = self.isrealin;
         % eigsopts.fail = 'keep';

         if inexact
             self.eigWork.k = min(k * 2, 20);
             [V, d, self.eigWork] = eig_inexact(Afun, self.eigWork, self.n);
         else
             [V,d,~] = eigs(Afun,self.n,k,sigma,eigsopts);
             d = diag(d);
         end
         % d = real(d);
         [d,permutation] = sort(real(d),'descend');
         V = V(:,permutation);
         f = d(1);
         if nargin >= 2
            g = self.forward(V(:,1));
         end
         self.nGDObjective = self.nGDObjective + 1;
      end % gdobjective

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [f,g,r] = pfdobjective(self,rhs,x)
%          r = self.forward(x)-rhs;
%          f = dot(r(:)/2,r(:)/2);
%          if nargout >= 2
%             g = self.adjoint(r,x);
%             if self.isrealin
%                g = real(g);
%             end
%             if nargout >= 3
%                r = -r;
%             end
%          end
         if (nargin <= 2 || isempty(x))
            if nargout == 1
               f = @(X)self.pfdobjective(rhs,X);
            else
               error('[hop.pl.pfdobjective] Invalid call!');
            end
         else
            x  = reshape(x,[self.n1,self.n2,1,numel(x)/self.n]);
            Ax = fft2(bsxfun(@times,self.masks,x));
            r  = rhs-dot(Ax,Ax,4)/self.n;
            f  = dot(r(:),r(:))/4;
            if nargout >= 2
               g  = -reshape(sum(bsxfun(@times,conj(self.masks),ifft2(bsxfun(@times,r,Ax))),3),[self.n1,self.n2,numel(x)/self.n]);
               if self.isrealin
                  g = real(g);
               end
            end
            util.addToNFFT2 (numel(Ax)/self.n);
            util.addToNIFFT2(numel(Ax)/self.n);
            self.nPFDObjective = self.nPFDObjective + 1;
         end
      end % pfdobjective

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [f,g,r] = dfpobjective(self,xlhs,xrhs,y)
         %DFPOBJECTIVE  Dual from primal.
%          Actyxlhs = self.adjoint(y,xlhs);
%          r = Actyxlhs(:)-xrhs(:);
%          f = dot(r(:),r(:))/2;
%          if nargout >= 2
%             g = self.forward(r,xlhs);
%             if self.isrealin
%                g = real(g);
%             end
%             if nargout >= 3
%                r = -r;
%             end
%          end
         xlhs = reshape(xlhs,[self.n1,self.n2,1,numel(xlhs)/self.n]);
         V    = fft2(bsxfun(@times,self.masks,xlhs));
         util.addToNFFT2(numel(V)/self.n);
         if (nargin <= 3 || isempty(y))
            if nargout == 1
               f = @(Y)self.dfpobjectiveV(V,xrhs,Y);
            else
               error('[hop.pl.dfpobjective] Invalid call!');
            end
         else
            switch nargout
               case 1
                  f = self.dfpobjectiveV(V,xrhs,y);
               case 2
                  [f,g] = self.dfpobjectiveV(V,xrhs,y);
               case 3
                  [f,g,r] = self.dfpobjectiveV(V,xrhs,y);
            end
         end
      end % dfpobjective

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function sz = size(self,dim)
         sz = [self.m,self.n,self.n];
         if nargin == 2
            if dim > 3
               sz = 1;
            elseif dim < 1
               error('pl:size',['Error using <strong>size</strong>\n'  ...
                     'Dimension argument must be a positive integer\n' ...
                     'scalar within indexing range.']);
            else
               sz = sz(dim);
            end
         end
      end % size

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   end % methods

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   methods ( Access = private )

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [f,g,r] = dfpobjectiveV(self,V,xrhs,y)
         xrhs = reshape(xrhs,[self.n1,self.n2,1      ,numel(xrhs)/self.n]);
         y    = reshape(y   ,[self.m1,self.m2,self.m3,1                 ]);
         if numel(xrhs) == self.n
            r = dot(self.masks,ifft2(y.*V),3)-xrhs;
         else
            r = sum(bsxfun(@times,conj(self.masks),ifft2(bsxfun(@times,y,V))),3)-xrhs;
         end
         f = dot(r(:),r(:))/2;
         util.addToNIFFT2(self.m3*numel(r)/self.n);
         if nargout >= 2
            Z = fft2(bsxfun(@times,self.masks,r));
            g = dot(V,Z,4)/self.n;
            util.addToNFFT2(numel(Z)/self.n);
            % if self.isrealin
               g = real(g);
            % end
            if nargout >= 3
               r = -r;
            end
         end
         self.nDFPObjective = self.nDFPObjective + 1;
      end % dfpobjectiveV

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   end % methods ( Access = private )

end % classdef pl < hop.hop
