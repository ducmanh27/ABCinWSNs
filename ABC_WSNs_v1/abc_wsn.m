clc;
clear;
close all;
N = 20 ; %Number of sensor nodes

Area = [100 100] %Area [x1 x2 .. xn y1 .. yn]

nVar = N * 2;             % Variables' dimensions

VarSize = [1 nVar];   % Decision Variables Matrix Size

VarMin = 0;         % Decision Variables Lower Bound

VarMax = 100;         % Decision Variables Upper Bound

tRange = 10 %Transmission range = 10 meter
%% ABC Settings

MaxIt = 500;              % Maximum Number of Iterations

nPop = 50;               % Population Size (Colony Size)

nOnlooker = nPop;         % Number of Onlooker Bees

Limit = nPop * nVar; % Abandonment Limit Parameter (Trial Limit)

a = 1;                    % Acceleration Coefficient Upper Bound

Ran = 10;
%% Random Initialization and Waggle dance
% Empty Bee Structure
empty_bee.Position = [];
empty_bee.Cost = [];

% Initialize Population Array
pop = repmat(empty_bee, nPop, 1);

% Initialize Best Solution Ever Found
BestSol.Cost = inf;

% Create Initial Population
for k = 1:nPop
    pos = zeros(N,2);
    pos (1,:)= [50 50];
    for i=2:N
        Rcom= Ran*rand;
        j=randi(i-1,1);
        pos (i,1)= pos(j,1)+2*Rcom*rand-Rcom;
        pos (i,2)= sqrt(Rcom^2-(pos(i,1)-pos(j,1))^2)+pos(j,2);
    end
    pop(k).Position =  reshape(pos, 1, []); %unifrnd(VarMin, VarMax, VarSize);
     CostFunction1 = Sphere(pop(k).Position, tRange, Area);
     pop(k).Cost = CostFunction1;
     if pop(k).Cost <= BestSol.Cost
         BestSol = pop(k);
     end
end 
% Abandonment Counter
 C = zeros(nPop, 1);
 
 % Array to Hold Best Cost Values
 BestCost = zeros(MaxIt, 1);
 BestFit = zeros(MaxIt,1);
 %% Local search
for it = 1:MaxIt
    for i = 1:nPop
        while (1)
            % Choose k randomly, not equal to i
            K = [1:i-1 i+1:nPop];
            k = K(randi([1 numel(K)])); %Randomly neibourhood bee
        
            % Define Acceleration Coeff.
            phi = a*unifrnd(-1, +1);
            w1 = randi([1 N*2]);
            % New Bee Position
            if w1>N
                w1_new=w1-N;
            else
                w1_new = w1+N;
            end
            newbee.Position = pop(i).Position;
            newbee.Position(w1) = pop(i).Position(w1)+phi.*(pop(i).Position(w1)-pop(k).Position(w1));
            newbee.Position(w1_new) = pop(i).Position(w1_new)+phi.*(pop(i).Position(w1_new)-pop(k).Position(w1_new));
            % Apply Bounds
            newbee.Position = max(newbee.Position, VarMin);
            newbee.Position = min(newbee.Position, VarMax);
            check1 = Connectivity(newbee.Position, N, Ran);
            if check1 == 0
                newbee.Position = pop(i).Position;
            else
                break;
            end
        end
    end
    
    % Calculate Fitness Values and Selection Probabilities
    F = zeros(nPop, 1);
    MeanCost = mean([pop.Cost]);
    for i = 1:nPop
        %-------------------------------------------------------------------
        F(i) = exp(-pop(i).Cost/MeanCost); % Convert Cost to Fitness
        %-------------------------------------------------------------------
    end
    P = F/sum(F);
    
    % Onlooker Bees
    for m = 1:nOnlooker % run m times but no use m as a matrix counter
        
        % Select Source Site
        %-------------------------------------------------------------------
        i = RouletteWheelSelection(P); %not run i from 1 to pop but randomly select i using roullette/greedy
        %-------------------------------------------------------------------
        % Choose k randomly, not equal to i
        while (1)
            K = [1:i-1 i+1:nPop];
            k = K(randi([1 numel(K)]));
        
        % Define Acceleration Coeff.
            phi = a*unifrnd(-1, +1);
            w2 = randi([1 N*2]);
            if w2>N
                w2_new=w2-N;
            else
                w2_new = w2+N;
            end
        % New Bee Position
            newbee.Position = pop(i).Position; 
            newbee.Position(w2) = pop(i).Position(w2)+phi.*(pop(i).Position(w2)-pop(k).Position(w2));
            newbee.Position(w2_new) = pop(i).Position(w2_new)+phi.*(pop(i).Position(w2_new)-pop(k).Position(w2_new));
        % Apply Bounds
             newbee.Position = max(newbee.Position, VarMin);
             newbee.Position = min(newbee.Position, VarMax);
            check1 = Connectivity(newbee.Position, N, Ran);
            if check1 == 0
                newbee.Position = pop(i).Position;
            else
                break;
            end
        end

        % Evaluation
        CostFunction2 = Sphere(newbee.Position, tRange, Area);
        newbee.Cost = CostFunction2; 
        
        % Comparision
        if newbee.Cost <= pop(i).Cost
            pop(i) = newbee;
        else
            C(i) = C(i) + 1;
        end
        
    end
    %-----------------------------------------------------------------------
    %% Global Search
    % Scout Bees
    for i = 1:nPop
        if C(i) >= Limit
            pos = zeros(N,2);
            pos (1,:)= [50 50];
            for m=2:N
                Rcom= Ran*rand;
                j=randi(m-1,1);
                pos (m,1)= pos(j,1)+2*Rcom*rand-Rcom;
                pos (m,2)= sqrt(Rcom^2-(pos(m,1)-pos(j,1))^2)+pos(j,2);
            end
            pop(i).Position =  reshape(pos, 1, []); %unifrnd(VarMin, VarMax, VarSize);
            CostFunction3 = Sphere(pop(k).Position, tRange, Area);
            pop(k).Cost = CostFunction3;
            C(i) = 0;
        end
    end
    %-----------------------------------------------------------------------delete
    %exhausted source out of loop
    
    % Update Best Solution Ever Found
    for i = 1:nPop
        if pop(i).Cost <= BestSol.Cost
            BestSol = pop(i);
        end
    end
    % Store Best Cost Ever Found
    BestCost(it) = BestSol.Cost;
    BestPoisition = BestSol.Position;
    BestFit(it) = exp(-BestSol.Cost/mean([pop.Cost]));
    % Display Iteration Information
    disp(['Iteration ' num2str(it) ': Best ratio coverage = ' num2str(100/BestCost(it)) '%']);   
end
%for i = 1:nPop
    %disp([ 'x' : num2str(pop(i).Position)]);
%end
%% Results
%pop.Position=categorical(pop.Position);
figure;
%plot(BestCost, 'LineWidth', 2);
%hold on;
semilogy(BestCost, 'LineWidth', 2);
%geobubble(pop.Position(1,:),pop.Position(2,:));
xlabel('Iteration');
ylabel('Best Fit');
 grid on;
 finalPos = reshape(BestSol.Position,[numel(BestSol.Position)/2,2]);
 
 figure
 plot(finalPos(:,1),finalPos(:,2),'o','color','r');
 hold on
 for ii=1:N                 % plot the circular transmission range
     [finalcircle.x(ii,:),finalcircle.y(ii,:)]=circle(finalPos(ii,1),finalPos(ii,2),tRange);
     fill(finalcircle.x(ii,:),finalcircle.y(ii,:),[0.25,0.25,0.25]);
     alpha 0.3
     hold on
 end
 axis on
 xlabel('x(m)')
 ylabel('y(m)')
 title('Optimized location of Nodes using Artificial bee colony algorithm')
 figure
 adj_matrix_final = zeros(N,N);
 for i=1:N
    for j=1:N
        if (((finalPos(i,1)-finalPos(j,1))^2+(finalPos(i,2)-finalPos(j,2))^2)<=(Ran)^2)
            adj_matrix_final(i,j) = 1;
        end
    end
 end
for i=1:N
    adj_matrix_final(i,i) = 0;
end
G= graph(adj_matrix_final);
plot(G);