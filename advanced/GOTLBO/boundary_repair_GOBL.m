function u = boundary_repair_GOBL(v,low,up,a,b)

[NP, D] = size(v);   
  
for i = 1:NP    
    for j = 1:D 
    if v(i,j) > up(j) || v(i,j) < low(j)
        u(i,j) = a(j) + rand*(b(j)-a(j));
    else
        u(i,j) = v(i,j);
    end  
end
end   
 