function V = initial_velocity(N,D);
[a,b]= find_boundary(D);
%a = -1;
%b =  1;
for i=1:N 
V(:,i)=(b-a).*rand(1,D)'*0.1;
end

