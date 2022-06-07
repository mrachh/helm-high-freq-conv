
function submat= helm_circ_kern(zk,srcinfo,targinfo,eta)
%CHNK.HELM2D.KERN standard Helmholtz layer potential kernels in 2D
% 
% Syntax: submat = helm_circ_kern(zk,srcinfo,targingo,eta)
%
% Let x be targets and y be sources for these formulas, with
% n_x and n_y the corresponding unit normals at those points
% (if defined). Note that the normal information is obtained
% by taking the perpendicular to the provided tangential deriviative
% info and normalizing  
%  
% Kernels based on G(x,y) = i/4 H_0^{(1)}(zk |x-y|)
%
% D(x,y) = \nabla_{n_y} G(x,y)
% S(x,y) = G(x,y)
% submat = D + eta S
%
% Input:
%   zk - complex number, Helmholtz wave number
%   srcinfo - description of sources in ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targinfo - description of targets in ptinfo struct format,
%                if info not relevant (d/d2) it doesn't need to
%                be provided. sprime requires tangent info in
%                targinfo.d
%   type - string, determines kernel type
%                type == 'd', double layer kernel D
%                type == 's', single layer kernel S
%                type == 'sprime', normal derivative of single
%                      layer S'
%                type == 'c', combined layer kernel D + i eta S
%   varargin{1} - eta in the combined layer formula, otherwise
%                does nothing
%
% Output:
%   submat - the evaluation of the selected kernel for the
%            provided sources and targets. the number of
%            rows equals the number of targets and the
%            number of columns equals the number of sources  
%
% see also CHNK.HELM2D.GREEN
  
src = srcinfo.r;
targ = targinfo.r;

tsrc = atan2(src(1,:),src(2,:));
ttarg = atan2(targ(1,:),targ(2,:));

[~,ns] = size(src);
[~,nt] = size(targ);

tsrc_rep = repmat(tsrc,[nt,1]);
ttarg_rep = repmat(ttarg.',[1,ns]);

tdiff = tsrc_rep-ttarg_rep;
tdiff(tdiff>2*pi) = tdiff(tdiff>2*pi) - 2*pi;
tdiff(tdiff<-2*pi) = tdiff(tdiff<-2*pi) + 2*pi;

z = 2*zk*abs(sin(tdiff/2));
h0 = besselh(0,1,z);
h1 = besselh(1,1,z);

submat = -1j/4*zk*abs(sin(tdiff/2)).*h1 + 1j*eta*(1j/4*h0);
