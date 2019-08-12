clear all; clc;

fprintf('  Seed          Best Cost      Iteration      cpu time\n')
percent=0;
seed=3;
n = 15;% number of variable
Sbest = zeros(n,seed);
maxtimerun=500;
for loop=1:seed,
      fprintf('   %d',loop)
      start_time = cputime;
      use_time=0;
      find_time=0;
      iteration=0;
      Fmin = 1;
      clear x F Fbest Fb cbest C W V gbest
      %Initialization of PSO parameters
      wmax=0.9;
      wmin=0.4;
      itmax=30; %Maximum iteration number
      c1 = 1.4;% constant of acceletion 
      c2 = 1.4;
      for iter=1:itmax
           W(iter)=wmax-((wmax-wmin)/itmax)*iter; % Inertial weight factor
      end
      iter=1;
      %B = zeros(1,1,itmax+1)
%**********************************************************
%Initialization of positions of agents
% agents are initialized between -5,+5 randomly
      m = 100;% number of agents
      [B_low,B_up] = find_boundary(n);      
      x = generate_initial(B_low,B_up,m,n);
%Initialization of velocities of agents
      V = initial_velocity(m,n);
%**********************************************************
%Function to be minimized.
     for i=1:m;
          con = check_constraint(x(:,i,1),n);
          if con == 1,
              F(i,1,1)=find_cost(x(:,i,1));
          else
              F(i,1,1)= 50000;
          end
     end
%**********************************************************
      [C,I]=min(abs(F(:,1,1)));
      B(1,1,1)=C;% min value
      XX(1,1,1)=I; % position of min value
      gbest(:,1,1)=x(:,I,1);% x at min value
%********************************************************
%Matrix composed of gbest vector 
      for p=1:m
          G(:,p,1)=gbest(:,1,1);
      end
      Fbest(1,1,1)=find_cost(G(:,1,1));
      for i=1:m;
           pbest(:,i,1)=x(:,i,1);
      end
      V(:,:,2)=W(1)*V(:,:,1)+c1*rand*(pbest(:,:,1)-x(:,:,1))+c2*rand*(G(:,:,1)-x(:,:,1));
      x(:,:,2)=x(:,:,1)+V(:,:,2);
      Fb(1,1,1)=find_cost(gbest(:,1,1));
      gx = 0;
%******************************************************
      while (iter < itmax) 
            iter = iter+1;
            iteration=iteration+1;
            j=iter;
                    % Calculation of new positions
                for i=1:m;
                    con = check_constraint(x(:,i,j),n);
                    if con==1,
                          F(i,1,j)=find_cost(x(:,i,j));
                    else
                          F(i,1,j)= 50000;
                    end
                end                
                         [C,I]=min(abs(F(:,:,j)));
                    %best position 
                    B(1,1,j)=C;            %lowest cost in round 
                    gbest(:,1,j)=x(:,I,j); %best position of lowest cost in round 
                    %
                            Fb(1,1,j)=find_cost(gbest(:,1,j));
                            [C,I]=min(Fb(1,1,:));
                            iterplot(iteration)=iteration;
                            BestCplot(iteration)=C;
                                  if Fb(1,1,j)<=C
                                       gbest(:,1,j)=gbest(:,1,j);
                                  else
                                       gbest(:,1,j)=gbest(:,1,I);
                                  end
                                       cbest = C;
      if mod(iter,1)==0
           sprintf('%d move, Best cost = %g',iter,C);
      end
    %Matrix composed of gbest vector 
      for p=1:m
           for r=1:n  
               G(:,p,j)=gbest(:,1,j);
           end
      end
                Fbest(1,1,j)=find_cost(G(:,1,j));
      for i=1:m;
                [C,I]=min(F(i,1,:));
          if F(i,1,j)<=C
              pbest(:,i,j)=x(:,i,j);
          else
              pbest(:,i,j)=x(:,i,I);
          end
      end
                       V(:,:,j+1)=W(j)*V(:,:,j)+c1*rand*(pbest(:,:,j)-x(:,:,j))+c2*rand*(G(:,:,j)-x(:,:,j));
            for t=1:m,
                    x(:,t,j+1)=x(:,t,j)+V(:,t,j+1);
                    con=check_constraint(x(:,t,j+1),n);
                    if con == 0
                            x(:,t,j+1)=x(:,t,j);
                    end
            end
            
                error = abs(cbest-Fmin);
                    if gx==0,
                            if error <= 1e-6,
                                stop_time = cputime;
                                find_time = stop_time - start_time;
                                percent=percent+1;
                                gx = 1;
                            end
                    end
                        check_time = cputime;
                        count_time = check_time - start_time;
                if count_time >= maxtimerun,
                       iter = itmax; 
                end
               
                
      end % while loop
                    end_time = cputime;
                    use_end_time = end_time - start_time;
                    a=size(gbest);                    
                    s_best = gbest(:,1,a(3));
                    Sbest(:,loop)=s_best;
                    best_c = find_cost(s_best);
                    fprintf('            %g',best_c)
                    error = abs(best_c-Fmin);
                if error <= 1e-6,
                    use_time = find_time;
                else
                    use_time=use_end_time;
                end
                        %use_time=use_end_time;
                        fprintf('          %g',iteration)
                        fprintf('          %g\n',use_time)
                        %Record
                        cost_r(loop)=best_c;
                        iter_r(loop)=iteration;
                        usetime_r(loop)=use_time;
end % for loop=1:seed
Percent=percent*100/seed;
[min_cost,min_p]=min(cost_r);
mean_cost=mean(cost_r);
max_cost=max(cost_r);
std_cost=std(cost_r);
cof_cost=std_cost/mean_cost;
mean_iter=mean(iter_r);
mean_cpu=mean(usetime_r);

fprintf('Cost           min =%g , average =%g,  max = %g ,  cof.variation = %g\n',min_cost,mean_cost,max_cost,cof_cost)
fprintf('Iteration      average =%g \n',mean_iter)
%fprintf('Search space   average =%g \n',mean_ss)
fprintf('CPU time       average =%g \n',mean_cpu)
fprintf('Percent get optimal    =%g \n',Percent)
[Minimum_Cost,Pg,Trm_loss,MVA_line]=find_cost_final(Sbest(:,min_p));
Result_caseieee30bus=[Pg(1);Sbest(:,min_p)];

fprintf('Variable    Value \n')
fprintf('PG1(MW) %7.3f\n', Result_caseieee30bus(1))
fprintf('PG2(MW) %7.3f\n', Result_caseieee30bus(2))
fprintf('PG5(MW) %7.3f\n', Result_caseieee30bus(3))
fprintf('PG8(MW) %7.3f\n', Result_caseieee30bus(4))
fprintf('PG11(MW) %7.3f\n', Result_caseieee30bus(5))
fprintf('PG13(MW) %7.3f\n', Result_caseieee30bus(6))
fprintf('VG1(pu) %7.3f\n', Result_caseieee30bus(7))
fprintf('VG2(pu) %7.3f\n', Result_caseieee30bus(8))
fprintf('VG5(pu) %7.3f\n', Result_caseieee30bus(9))
fprintf('VG8(pu) %7.3f\n', Result_caseieee30bus(10))
fprintf('VG11(pu) %7.3f\n', Result_caseieee30bus(11))
fprintf('VG13(pu) %7.3f\n', Result_caseieee30bus(12))
fprintf('T11 %7.3f\n', Result_caseieee30bus(13))
fprintf('T12 %7.3f\n', Result_caseieee30bus(14))
fprintf('T15 %7.3f\n', Result_caseieee30bus(15))
fprintf('T36 %7.3f\n', Result_caseieee30bus(15))
fprintf('Fuel cost($/h) %7.3f\n', min_cost)
fprintf('Transmissionloss(MW) %7.3f\n', Trm_loss)
fprintf('CPU time(s) %7.3f\n', usetime_r(:,min_p))
plot(iterplot,BestCplot)
