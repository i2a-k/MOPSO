clc;
clear;
close all;

pareto_front_ZDT = xlsread('ZDT1.xlsx');

%Number of Initial Population & Dimensions
pop_size = 200;
dim = 2;

%Number of Objectives
k = 2;

from = 0;
to = 1;

%Maximum Number of Answers for Final Pareto Front
Pareto_Front = zeros(pop_size,dim+k);

%Iteration Condition
max_iter = 100;

%Results for n Times Execution
num_of_result = 5;
%Total Results
total_time = zeros(num_of_result,1);
GD_A = zeros(num_of_result,1);
Spread_A = zeros(num_of_result,1);
total_Pareto_Front = zeros(pop_size,dim+k,num_of_result);

for n=1:num_of_result
    tic;
    
    %Initial Global Best
    g_best = zeros(1, dim);
    g_best_fitness = inf;
    
    %Initial Personal Bests
    p_best = zeros(pop_size,dim);
    p_best(:,:) = 100;
    
    %Initial Population
    X = unifrnd(from,to,[pop_size dim]);
    
    %Initial Velocity of each Particle
    %V = X;
    V = zeros(pop_size,dim);
    
    F_result = zeros(1, pop_size);
    
    c_max = 2.5;
    c_min = 0.5;
    
    for j = 1:max_iter
        %Calculate Fitness of X and Update Pbest and Gbest
        for i = 1:pop_size
            F_result(1,i) = AOF(X(i,:),k);
            if (F_result(1,i) <= AOF(p_best(i,:),k))
                p_best(i,:) = X(i,:);
            end
            
            if (F_result(1,i) <= g_best_fitness)
                g_best_fitness = F_result(1,i);
                g_best(1,:) = X(i,:);
            end
        end
        
        r1 = rand;
        r2 = rand;
        
        %Cognitive Coefficient
        c1 = (c_min - c_max) * (j/max_iter) + c_max;
        %Global Coefficient
        c2 = (c_max - c_min) * (j/max_iter) + c_max;
        %Inertia Coefficient for Standard PSO
        W = c1*r1 + c2*r2;
        
        X_new = zeros(pop_size,dim);
        
        for l=1:pop_size
            %Equation of Velocity (Update Velocity of each Particle)
            V(l,:) = (W * V(l,:)) + (c1*r1*(p_best(l,:) - X(l,:)))...
                + (c2*r2*(g_best(1,:) - X(l,:) )); %Standard PSO
            
            %Control the Domain of the new Velocity
            for p=1:dim
                r_v = rand;
                if (V(l,p) < from)
                    V(l,p) = from + r_v;
                elseif (V(l,p) > to)
                    V(l,p) = to - r_v;
                end
            end
            
            %Equation of new Position (Update Position of each Particle)
            X_new(l,:) = X(l,:) + V(l,:);
            
            %Control the Domain of the new Positions
            for t=1:dim
                r_x = rand;
                if (X_new(l,t) < from)
                    X_new(l,t) = from + r_x;
                elseif (X_new(l,t) > to)
                    X_new(l,t) = to - r_x;
                end
            end
        end
        
        %Combine X and X_new and Create New Swarm 
        NX = cat(1,X,X_new); %(NX = New Swarm)
        NX_F = zeros(size(NX,1),dim+k); %(NX_F = Fitness of New Swarm)
        NX_F(:,1:dim) = NX(:,1:dim);
        
        %Calculate Fitness of each Particle on ZDT
        for l=1:size(NX,1)
            NX_F(l,dim+1:end) = ZDT1(NX(l,:));
        end
        
        %Non-Dominating Sorting
        %if p Dominate q => add q to SP
        %if q Dominate p => (np=np+1)
        %if np=0 => add p to FI 
        SP = [];
        np = 0;
        % FI = [];
        FI = zeros(size(NX_F,1),dim+k);
        FI_Counter = 1;
        for p = 1:size(NX_F,1)-1
            for q = p+1:size(NX_F,1)
                if Dominates(NX_F(p,dim+1:end),NX_F(q,dim+1:end))
                    if( ~all(ismember(NX_F(q,:),SP)) )
                        SP = cat(1,SP,NX_F(q,:));
                    end
                elseif Dominates(NX_F(q,dim+1:end),NX_F(p,dim+1:end))
                    np = np + 1;
                end
                
                if np == 0
                    if( ~all(ismember(NX_F(p,:),FI)) )
                        %FI = cat(1,FI,NX_F(p,:));
                        FI(FI_Counter,:) = NX_F(p,:);
                        FI_Counter = FI_Counter + 1;
                    end
                end
            end
        end
        i = 1;
        while ( ~all(FI(i,:) == 0) )
            H = [];
            nq = size(SP,1);
            for p = 1:size(FI,1)
                for q = 1:size(SP,1)
                    nq = nq - 1;
                    %if nq=0 add q to H
                    if nq == 0
                        H = cat(1,H,SP(q,1:dim+k));
                    end
                end
            end
            i = i + 1;
            FI(i,:) = H;
            FI_Counter = FI_Counter + 1;
            i = i + 1;
        end
        
        %Select New Population for Next Iteration
        %Picking up FI Answers
        if( (FI_Counter - 1) > 0 )
            X(1:FI_Counter-1,:) = FI(1:FI_Counter-1,1:dim);
            Pareto_Front(1:FI_Counter-1,:) = FI(1:FI_Counter-1,:);
        end
        %Picking up SP Answers
        SP = flipud(SP);
        X(FI_Counter:end,:) = SP(1:pop_size-FI_Counter+1,1:dim);
        Pareto_Front(FI_Counter:end,:) = SP(1:pop_size-FI_Counter+1,:);
        
        figure(n)
        plot(Pareto_Front(:,dim+(k-1)),Pareto_Front(:,dim+k),'ko')
        hold on
        plot(pareto_front_ZDT(:,1),pareto_front_ZDT(:,2),'r*')
        xlabel('f1');
        ylabel('f2');
        grid on;
        hold off;
        pause(0.001);
        
        disp( strcat( 'Iteration: ', num2str(j) ) )
    end
    
    total_time(n,1) = toc;
    total_Pareto_Front(:,:,n) = Pareto_Front(:,:);
    [di, dm] = Di_Dm(Pareto_Front(:,dim+1:end),pareto_front_ZDT);
    [GD_A(n,1), Spread_A(n,1)] = GD_SA(di,size(Pareto_Front,1),dm,Pareto_Front(:,dim+1:end));
    
end

mean_GD = mean(GD_A(:,1));
mean_SA = mean(Spread_A(:,1));
std_GD = std(GD_A(:,1));
std_SA = std(Spread_A(:,1));
mean_time = mean(total_time(:,1));

disp(strcat('mean GD: ', num2str(mean_GD)));
disp(strcat('mean SA: ', num2str(mean_SA)));
disp(strcat('std generational distance: ', num2str(std_GD)));
disp(strcat('std spread assessment: ', num2str(std_SA)));
disp(strcat('mean time: ', num2str(mean_time)));
