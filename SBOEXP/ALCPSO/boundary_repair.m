function u = boundary_repair(v,low,up,str)

[NP, D] = size(v);   
u = v; 
llow = repmat(low,NP,1);
uup  = repmat(up,NP,1);

if strcmp(str,'absorb')
    index1 = (u>uup);
    u(index1) = uup(index1); 
    index2 = (u<llow);
    u(index2) = llow(index2); 
end
   

if strcmp(str,'random')
    index1 = (u>uup)|(u<llow);
    rr = rand(size(u));
    u(index1) = llow(index1) + rr(index1).*(uup(index1)-llow(index1));  
end


if strcmp(str,'reflect')
    index1 = (u>uup);
    u(index1) = max( 2*uup(index1)-u(index1), llow(index1) );
    index2 = (u<llow);
    u(index2) = min( 2*llow(index2)-u(index2), uup(index2) );
end
