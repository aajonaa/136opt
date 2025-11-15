% Harris' Hwak Algorithm modified by Jona 2023-10-27.
function [best_pos,convergence_curve]= OLHHO(N,Max_FEs,lb,ub,dim,fobj)
tic
disp('OLHHO is now tackling your problem');
% initialize the location and Energy of the rabbit
best_pos=zeros(1,dim);
bestFitness=inf;

%Initialize the locations of Harris' hawks
X=initialization(N,dim,ub,lb);

convergence_curve=[];

t=1; % Loop counter

FEs = 0;
while FEs < Max_FEs
    for i=1:size(X,1)
        % Check boundries
        FU=X(i,:)>ub;FL=X(i,:)<lb;X(i,:)=(X(i,:).*(~(FU+FL)))+ub.*FU+lb.*FL;
        % fitness of locations
        fitness=fobj(X(i,:));
        FEs = FEs + 1;
        % Update the location of Rabbit
        if fitness<bestFitness
            bestFitness=fitness;
            best_pos=X(i,:);
        end
    end
    
    E1=2*(1-(FEs/Max_FEs)); % factor to show the decreaing energy of rabbit
    % Update the location of Harris' hawks
    for i=1:size(X,1)
        E0=2*rand()-1; %-1<E0<1
        Escaping_Energy=E1*(E0);  % escaping energy of rabbit
        
        if abs(Escaping_Energy)>=1
            %% Exploration:
            % Harris' hawks perch randomly based on 2 strategy:
            
            q=rand();
            rand_Hawk_index = floor(N*rand()+1);
            X_rand = X(rand_Hawk_index, :);
            if q<0.5
                % perch based on other family members
                X(i,:)=X_rand-rand()*abs(X_rand-2*rand()*X(i,:));
            elseif q>=0.5
                % perch on a random tall tree (random site inside group's home range)
                X(i,:)=(best_pos(1,:)-mean(X))-rand()*((ub-lb)*rand+lb);
            end
            
        elseif abs(Escaping_Energy)<1
            %% Exploitation:
            % Attacking the rabbit using 4 strategies regarding the behavior of the rabbit
            
            %% phase 1: surprise pounce (seven kills)
            % surprise pounce (seven kills): multiple, short rapid dives by different hawks
            
            r=rand(); % probablity of each event
            
            if r>=0.5 && abs(Escaping_Energy)<0.5 % Hard besiege
                X(i,:)=(best_pos)-Escaping_Energy*abs(best_pos-X(i,:));
            end
            
            if r>=0.5 && abs(Escaping_Energy)>=0.5  % Soft besiege
                Jump_strength=2*(1-rand()); % random jump strength of the rabbit
                X(i,:)=(best_pos-X(i,:))-Escaping_Energy*abs(Jump_strength*best_pos-X(i,:));
            end
            
            %% phase 2: performing team rapid dives (leapfrog movements)
            if r<0.5 && abs(Escaping_Energy)>=0.5 % Soft besiege % rabbit try to escape by many zigzag deceptive motions
                
                Jump_strength=2*(1-rand());
                X1=best_pos-Escaping_Energy*abs(Jump_strength*best_pos-X(i,:));
                
                FEs = FEs + 1;
                if fobj(X1)<fobj(X(i,:)) % improved move?
                    X(i,:)=X1;
                else % hawks perform levy-based short rapid dives around the rabbit
                    X2=best_pos-Escaping_Energy*abs(Jump_strength*best_pos-X(i,:))+rand(1,dim).*Levy(dim);
                    
                    FEs = FEs + 1;
                    if (fobj(X2)<fobj(X(i,:))) % improved move?
                        X(i,:)=X2;
                    end
                end
            end
            
            if r<0.5 && abs(Escaping_Energy)<0.5 % Hard besiege % rabbit try to escape by many zigzag deceptive motions
                % hawks try to decrease their average location with the rabbit
                Jump_strength=2*(1-rand());
                X1=best_pos-Escaping_Energy*abs(Jump_strength*best_pos-mean(X));
                
                FEs = FEs + 1;
                if fobj(X1)<fobj(X(i,:)) % improved move?
                    X(i,:)=X1;
                else % Perform levy-based short rapid dives around the rabbit
                    X2=best_pos-Escaping_Energy*abs(Jump_strength*best_pos-mean(X))+rand(1,dim).*Levy(dim);
                    
                    FEs = FEs + 1;
                    if (fobj(X2)<fobj(X(i,:))) % improved move?
                        X(i,:)=X2;
                    end
                end
            end
        end
    end
    
    %%%%%%%%%%
    dynamic_ub=max(X);
    dynamic_lb=min(X);
     for i = 1:N  
        S(i,:) = rand() * (dynamic_ub+dynamic_lb) - X(i,:);
         for j = 1:dim
            if S(i,j) < lb || S(j) > ub
                S(i,j) = lb + rand() * (ub - lb);  
            end
         end 
     end 
    Positions1 = S;
    Positions2 = [X; Positions1];
    for i = 1:2*N
        fit_value(i) = fobj(Positions2(i, :));
    end
    [fit_value,index] = sort(fit_value);
    P = Positions2(index, :); 
    X = P(1:N, :); 
    %%%%%%%%%%
    
    convergence_curve(t)=bestFitness;
    t=t+1;
end
toc
end

function o=Levy(d)
beta=1.5;
sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
u=randn(1,d)*sigma;v=randn(1,d);step=u./abs(v).^(1/beta);
o=step;
end

