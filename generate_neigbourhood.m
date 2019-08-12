function x = generate_neigbourhood(N_size,xr,nep,n);
s = zeros(1,n);
[B_low,B_up]= find_boundary(n);
%X_up  = xr + N_size'.*rand(n,1);
%X_low = xr - N_size'.*rand(n,1);
for k = 1:nep,
     s = zeros(n,1);
     Xn=xr+N_size(:,1)*(2*rand(n,1)-ones(n,1));
    for i=1:n,
        if Xn(i) >= B_up(i),
            Xn(i) = B_up(i);
        end
        if Xn(i) <= B_low(i),
            Xn(i)=B_low(i);
        end
     s(i)= Xn(i);
     x(:,k)= s;
    end     
end
