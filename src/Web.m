function [w]=Web(r0,r1,dw,L0)
    if (L0/2) < (r0-r1)
    w_max = L0/2;
    else
    w_max = r0-r1;
    end
    w=0:dw:w_max;
    if w(length(w))~=w_max
    w(length(w)+1)=w_max;
    end
    
end