function con = check_constraint(p,n);
[B_low,B_up]=find_boundary(n);
index = zeros(n,1);
for i=1:n,
    if (p(i) > B_low(i))&( p(i)< B_up(i))
        index(i)=1;
    end
end
con = prod(index);

    

