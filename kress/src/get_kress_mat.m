function [mat] = get_kress_mat(kh,src,t,spars,eta)
%
% Get kress matrix with combined field representation
%
 
    if(nargin <= 3)
        spars.ifsplit = false; 
    end
    
    if(nargin <= 4)
        eta = -kh;
    end
    
    S  = slmat(kh,src,t,spars);
    D  = dlmat(kh,src,t,spars);

    %Combined layer
    %solving for the field 
    
    [~,n] = size(src);
    
    mat = 2*(eye(n)/2 + D + 1i*eta*S);
    
    
end