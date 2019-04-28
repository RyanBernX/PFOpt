classdef hop < handle

   properties
      isrealin
      nForward
      nAdjoint
      nGDObjective
      nPFDObjective
      nDFPObjective
   end % properties
   
   methods ( Abstract )
      out       = forward(self,x1,x2,scatter)
      out       = adjoint(self,y,x)
      [f,g,V,d] = gdobjective(self,y,k,opts,v)
      [f,g,r]   = pfdobjective(self,rhs,x)
      [f,g,r]   = dfpobjective(self,xlhs,xrhs,y)
      sz        = size(self,dim)
   end % methods ( Abstract )
   
   methods

      function self = hop(isrealin)
         if (nargin == 0) || isempty(isrealin) || ~islogical(isrealin)
            isrealin = false;
         end
         self.isrealin      = isrealin;
         self.nForward      = 0;
         self.nAdjoint      = 0;
         self.nGDObjective  = 0;
         self.nPFDObjective = 0;
         self.nDFPObjective = 0;
      end % hop (constructor)
      
   end % methods
   
end % classdef hop < handle
