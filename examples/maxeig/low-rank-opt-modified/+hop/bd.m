classdef bd < hop.hop

   properties
      b1
      b2
      fm
      m1
      m2
      m
      n1
      n2
      n
      i1
      i2
      eigWork
   end % properties

   methods

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function self = bd(b1,b2,isrealin,fm)
         assert((b1.m1 == b2.m1) && (b1.m2 == b2.m2));
         if nargin < 2
            isrealin = false;
         end
         self@hop.hop(isrealin);
         self.b1 = b1;
         self.b2 = b2;
         self.fm = fm;
         self.m1 = self.b1.m1;
         self.m2 = self.b1.m2;
         self.m  = self.m1*self.m2;
         self.n1 = size(self.b1,2);
         self.n2 = size(self.b2,2);
         self.n  = self.n1+self.n2;
         self.i1 =          1:self.n1 ;
         self.i2 = self.n1+(1:self.n2);
         self.eigWork = struct('force_exact', 5);
      end % bd (constructor)

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function out = forward(self,x1,x2,scatter)
         %FORWARD  Forward Blind Deconvolution operator.
         if (nargin < 4) || isempty(scatter)
            scatter = false;
         end
         x1 = reshape(x1,self.n,numel(x1)/self.n);
         if (nargin >= 3) && ~isempty(x2)
            x2 = reshape(x2,self.n,numel(x2)/self.n);
         else
            x2 = x1;
         end
         U1 =      fft2(self.b1.forward(     x1(self.i1,:) )) ;
         V2 = conj(fft2(self.b2.forward(conj(x2(self.i2,:)))));
         util.addToNFFT2(numel(U1)/self.m+numel(V2)/self.m);
         if scatter
            V2 = reshape(V2,[self.m1,self.m2,1,numel(V2)/self.m]);
         end
         out = bsxfun(@times,U1,conj(V2));
         if ~scatter
            out = sum(out,3);
         end
         if ~self.fm
            out = ifft2(out)/realsqrt(self.m);
            util.addToNIFFT2(numel(out)/self.m);
         else
            out = out/self.m;
         end
         self.nForward = self.nForward + 1;
      end % forward

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function out = adjoint(self,y,x)
         x = reshape(x,self.n,numel(x)/self.n);
         y = reshape(y,[self.m1,self.m2,1,numel(y)/self.m]);
         if ~self.fm
            Y = fft2(y)/realsqrt(self.m);
            util.addToNFFT2(numel(Y)/self.m);
         else
            Y = y;
         end
         out = reshape(.5*[     self.b1.adjoint(ifft2(bsxfun(@times,Y,conj(fft2(self.b2.forward(conj(x(self.i2,:)))))))) ; ...
                           conj(self.b2.adjoint(ifft2(bsxfun(@times,Y,conj(fft2(self.b1.forward(     x(self.i1,:) )))))))],...
                       [self.n,numel(x)/self.n,numel(y)/self.m]);
         util.addToNFFT2 (2*numel(x)/self.n);
         util.addToNIFFT2(2*numel(x)/self.n*numel(Y)/self.m);
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
         vec = @(x) x(:);
         % [begin] optimize FFT counts during eigensolve
         disableFM = ~self.fm;
         if ~self.fm
            y       = fft2(y)/realsqrt(self.m);
            self.fm = true;
            util.addToNFFT2(numel(y)/self.m);
         end
         % [ end ] optimize FFT counts during eigensolve
         if self.isrealin
            Afun  = @(x) reshape(real(vec(self.adjoint(y,x))), size(x));
            sigma = 'LA';
         else
            Afun  = @(x) reshape(vec(self.adjoint(y,x)), size(x));
            sigma = 'LR';
         end
         eigsopts.issym  = true;
         eigsopts.isreal = self.isrealin;

         if inexact
             self.eigWork.k = k * 2;
             [V, d, self.eigWork] = eig_inexact(Afun, self.eigWork, self.n);
         else
             [V,d] = eigs(Afun,self.n,k,sigma,eigsopts);
             d = diag(d);
         end
         % [begin] optimize FFT counts during eigensolve
         if disableFM
            self.fm = false;
         end
         % [ end ] optimize FFT counts during eigensolve
         [d,permutation] = sort(real(d),'descend');
         V = V(:,permutation);
         f = d(1);
         if nargin >= 2
            g = self.forward(V(:,1),V(:,1));
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
         if ~self.fm
            rhs = fft2(rhs)/realsqrt(self.m);
            util.addToNFFT2(numel(rhs)/self.m);
         end
         if (nargin <= 2 || isempty(x))
            if nargout == 1
               f = @(X)self.pfdobjectiveRHS(rhs,X);
            else
               error('[hop.bd.pfdobjective] Invalid call!');
            end
         else
            switch nargout
               case 1
                  f = self.pfdobjectiveRHS(rhs,x);
               case 2
                  [f,g] = self.pfdobjectiveRHS(rhs,x);
               case 3
                  [f,g,r] = self.pfdobjectiveRHS(rhs,x);
            end
         end
      end % pfdobjective

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      function [f,g,r] = dfpobjective(self,xlhs,xrhs,y)
%          Actyxlhs = self.adjoint(y,xlhs);
%          r = Actyxlhs(:)-xrhs(:);
%          f = dot(r(:),r(:))/2;
%          if nargout >= 2
%             g = .5*(self.forward(r,xlhs)+self.forward(xlhs,r));
%             if nargout >= 3
%                r = -r;
%             end
%          end
         xlhs = reshape(xlhs,self.n,numel(xlhs)/self.n);
         Ulhs =      fft2(self.b1.forward(     xlhs(self.i1,:) )) ;
         Vlhs = conj(fft2(self.b2.forward(conj(xlhs(self.i2,:)))));
         util.addToNFFT2(2*numel(xlhs)/self.n)
         if (nargin <= 3 || isempty(y))
            if nargout == 1
               f = @(Y)self.dfpobjectiveUV(Ulhs,Vlhs,xrhs,Y);
            else
               error('[hop.bd.dfpobjective] Invalid call!');
            end
         else
            switch nargout
               case 1
                  f = self.dfpobjectiveUV(Ulhs,Vlhs,xrhs,y);
               case 2
                  [f,g] = self.dfpobjectiveUV(Ulhs,Vlhs,xrhs,y);
               case 3
                  [f,g,r] = self.dfpobjectiveUV(Ulhs,Vlhs,xrhs,y);
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
               error('bd:size',['Error using <strong>size</strong>\n'  ...
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

      function [f,g,r] = pfdobjectiveRHS(self,RHS,x)
         x = reshape(x,self.n,numel(x)/self.n);
         U =      fft2(self.b1.forward(     x(self.i1,:) )) ;
         V = conj(fft2(self.b2.forward(conj(x(self.i2,:)))));
         r = dot(V,U,3)/self.m-RHS;
         f = dot(r(:),r(:))/4;
         if nargout >= 2
            g = [     self.b1.adjoint(ifft2(bsxfun(@times,r,     V ) )) ;...
                 conj(self.b2.adjoint(ifft2(bsxfun(@times,r,conj(U)))))]/2;
            if self.isrealin
               g = real(g);
            end
            if nargout >= 3
               r = -r;
               if ~self.fm
                  r = ifft2(r)*realsqrt(self.m);
                  util.addToNIFFT2(numel(r)/self.m);
               end
            end
            util.addToNIFFT2(2*numel(x)/self.n);
         end
         util.addToNFFT2(2*numel(x)/self.n);
         self.nPFDObjective = self.nPFDObjective + 1;
      end

      function [f,g,r] = dfpobjectiveUV(self,Ulhs,Vlhs,xrhs,y)
         xrhs = reshape(xrhs,self.n ,numel(xrhs)/self.n);
         y    = reshape(y   ,self.m1,self.m2           );
         if ~self.fm
            Y = fft2(y)/realsqrt(self.m);
            util.addToNFFT2(numel(Y)/self.m);
         else
            Y = y;
         end
         r = [     self.b1.adjoint(ifft2(bsxfun(@times,Y,     Vlhs ) )) ;...
              conj(self.b2.adjoint(ifft2(bsxfun(@times,Y,conj(Ulhs)))))]/2-xrhs;
         f = dot(r(:),r(:))/2;
         if nargout >= 2
            Ur =      fft2(self.b1.forward(     r(self.i1,:) )) ;
            Vr = conj(fft2(self.b2.forward(conj(r(self.i2,:)))));
            g  = (dot(Vlhs,Ur,3)+dot(Vr,Ulhs,3))/2;
            if ~self.fm
               g = ifft2(g)/realsqrt(self.m);
               util.addToNFFT2(numel(g)/self.m);
            else
               g = g/self.m;
            end
            if nargout >= 3
               r = -r;
            end
            util.addToNIFFT2(numel(Ur)/self.m+numel(Vr)/self.m);
         end
         util.addToNIFFT2((numel(Ulhs)/self.m+numel(Vlhs)/self.m)*(numel(Y)/self.m));
         self.nDFPObjective = self.nDFPObjective + 1;
      end % dfpobjectiveUV

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   end % methods ( Access = private )

end % classdef bd < hop.hop
