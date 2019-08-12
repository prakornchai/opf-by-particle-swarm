function x = generate_initial(B_low,B_up,m,n);
s = zeros(1,n);
for k = 1:m,
    con=0;
    while con==0
         for i=1:n
              %s(i)= B_low(i) + ((B_up(i)-B_low(i))*k/m);
              %s(i)= B_low(i) + ((B_up(i)-B_low(i))*rand(1)*k/m);
              %s(i)= B_low(i) + ((B_up(i)-B_low(i))*k/m*(2*rand(1)-1));
              s(i)= B_low(i) + ((B_up(i)-B_low(i))*rand(1));
         end
         con=check_constraint(s,n);
    end
    x(:,k)= s;
end
