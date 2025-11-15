function [best_pos,Convergence_curve]=PO(N,Max_FEs,lb,ub,dim,fobj)
    parties = 6;
    areas = 5;
    FEs = 0;
    tic
    % initialize position vector and score for the leader
    best_pos=zeros(1,dim);
    bestFitness=inf; %change this to -inf for maximization problems
    %Initialize the positions of search agents
    X=initialization(N,dim,ub,lb);
    auxPositions = X;
    prevPositions = X;
    Convergence_curve=[];
    AllFitness=inf * ones(N, 1);
    %Running phases for initializations
    % Election;   %Run election phase

    for i=1:size(X,1)
        lambda = 1 - FEs/Max_FEs;
        % Return back the search agents that go beyond the boundaries of the search space
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

        %Calculate objective function for each search agent
        AllFitness(i,1)=fobj(X(i,:));
        FEs = FEs + 1;

        %Update the leader
        if AllFitness(i,1)<bestFitness % Change this to > for maximization problem
            bestFitness=AllFitness(i,1); 
            best_pos=X(i,:);
        end        
    end

    auxFitness = AllFitness;
    prevFitness = AllFitness;
    % GovernmentFormation;

    aWinnerInd=zeros(areas,1);   %Indices of area winners in x
    aWinners = zeros(areas,dim); %Area winners are stored separately
    for a = 1:areas
        [aWinnerFitness,aWinnerParty]=min(AllFitness(a:areas:N));
        aWinnerInd(a,1) = (aWinnerParty-1) * areas + a;
        aWinners(a,:) = X(aWinnerInd(a,1),:);
    end    
    %Finding party leaders
    pLeaderInd=zeros(parties,1);    %Indices of party leaders in x
    pLeaders = zeros(parties,dim);  %X of party leaders in x
    for p = 1:parties
        pStIndex = (p-1) * areas + 1;
        pEndIndex = pStIndex + areas - 1;
        [partyLeaderFitness,leadIndex]=min(AllFitness(pStIndex:pEndIndex)); 
        pLeaderInd(p,1) = (pStIndex - 1) + leadIndex; %Indexof party leader
        pLeaders(p,:) = X(pLeaderInd(p,1),:);
    end

    t=0;% Loop counter
    while t<Max_FEs
        prevFitness = auxFitness;
        prevPositions = auxPositions;
        auxFitness = AllFitness;
        auxPositions = X;
    %     ElectionCampaign;   

        for whichMethod = 1:2
            for a = 1:areas
                for p = 1:parties
                    i = (p-1)*areas + a; %index of member

                    for j=1:dim
                        if whichMethod == 1         %position-updating w.r.t party leader
                            center = pLeaders(p,j);
                        elseif whichMethod == 2     %position-updating w.r.t area winner
                            center = aWinners(a,j);
                        end

                        %Cases of Eq. 9 in paper
                        if prevFitness(i) >= AllFitness(i) 
                            if (prevPositions(i,j) <= X(i,j) && X(i,j) <= center) ...
                                    || (prevPositions(i,j) >= X(i,j) && X(i,j) >= center)

                                radius = center - X(i,j); 
                                X(i,j) = center + rand() * radius;
                            elseif (prevPositions(i,j) <= X(i,j) && X(i,j) >= center && center >= prevPositions(i,j)) ...
                                    || (prevPositions(i,j) >= X(i,j) && X(i,j) <= center && center <= prevPositions(i,j))

                                radius = abs(X(i,j) - center);
                                X(i,j) = center + (2*rand()-1) * radius;
                            elseif (prevPositions(i,j) <= X(i,j) && X(i,j) >= center && center <= prevPositions(i,j)) ...
                                    || (prevPositions(i,j) >= X(i,j) && X(i,j) <= center && center >= prevPositions(i,j))

                                radius = abs(prevPositions(i,j) - center);
                                X(i,j) = center + (2*rand()-1) * radius;
                            end

                        %Cases of Eq. 10 in paper
                        elseif prevFitness(i) < AllFitness(i) 
                            if (prevPositions(i,j) <= X(i,j) && X(i,j) <= center) ...
                                    || (prevPositions(i,j) >= X(i,j) && X(i,j) >= center)


                                radius = abs(X(i,j) - center);
                                X(i,j) = center + (2*rand()-1) * radius;
                            elseif (prevPositions(i,j) <= X(i,j) && X(i,j) >= center && center >= prevPositions(i,j)) ...
                                    || (prevPositions(i,j) >= X(i,j) && X(i,j) <= center && center <= prevPositions(i,j))

                                radius = X(i,j) - prevPositions(i,j);
                                X(i,j) = prevPositions(i,j) + rand() * radius;
                            elseif (prevPositions(i,j) <= X(i,j) && X(i,j) >= center && center <= prevPositions(i,j)) ...
                                    || (prevPositions(i,j) >= X(i,j) && X(i,j) <= center && center >= prevPositions(i,j))

                                center2 = prevPositions(i,j);
                                radius = abs(center - center2);
                                X(i,j) = center + (2*rand()-1) * radius;
                            end
                        end

                    end
                end
            end
        end

        %     PartySwitching;

        psr = (1-t*((1)/Max_FEs)) * lambda;
        for p=1:parties
            for a=1:areas
                fromPInd = (p-1)*areas + a;
                if rand() < psr
                    %Selecting a party other than current where want to send the
                    %member
                    toParty = randi(parties);
                    while(toParty == p)
                        toParty = randi(parties);
                    end

                    %Deciding member in TO party
                     toPStInd = (toParty-1) * areas + 1;
                     toPEndIndex = toPStInd + areas - 1;
                     [~,toPLeastFit] = max(AllFitness(toPStInd:toPEndIndex));
                     toPInd = toPStInd + toPLeastFit-1;


                    %Deciding what to do with member in FROM party and switching
                    fromPInd = (p-1)*areas + a;
                    temp = X(toPInd,:);
                    X(toPInd,:) = X(fromPInd);
                    X(fromPInd,:)=temp;

                    temp = AllFitness(toPInd);
                    AllFitness(toPInd) = AllFitness(fromPInd);
                    AllFitness(fromPInd) = temp;
                end
            end
        end


    %     Election;

        for i=1:size(X,1)
            % Return back the search agents that go beyond the boundaries of the search space
            Flag4ub=X(i,:)>ub;
            Flag4lb=X(i,:)<lb;
            X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;

            %Calculate objective function for each search agent
            AllFitness(i,1)=fobj(X(i,:));
            FEs = FEs + 1;

            %Update the leader
            if AllFitness(i,1)<bestFitness % Change this to > for maximization problem
                bestFitness=AllFitness(i,1); 
                best_pos=X(i,:);
            end        
        end

    %     GovernmentFormation;

        aWinnerInd=zeros(areas,1);   %Indices of area winners in x
        aWinners = zeros(areas,dim); %Area winners are stored separately
        for a = 1:areas
            [aWinnerFitness,aWinnerParty]=min(AllFitness(a:areas:N));
            aWinnerInd(a,1) = (aWinnerParty-1) * areas + a;
            aWinners(a,:) = X(aWinnerInd(a,1),:);
        end    
        %Finding party leaders
        pLeaderInd=zeros(parties,1);    %Indices of party leaders in x
        pLeaders = zeros(parties,dim);  %X of party leaders in x
        for p = 1:parties
            pStIndex = (p-1) * areas + 1;
            pEndIndex = pStIndex + areas - 1;
            [partyLeaderFitness,leadIndex]=min(AllFitness(pStIndex:pEndIndex)); 
            pLeaderInd(p,1) = (pStIndex - 1) + leadIndex; %Indexof party leader
            pLeaders(p,:) = X(pLeaderInd(p,1),:);
        end

    %     Parliamentarism;


        for a=1:areas
            newAWinner = aWinners(a,:);
            i = aWinnerInd(a);    

            toa = randi(areas);
            while(toa == a)
                toa = randi(areas);
            end
            toAWinner = aWinners(toa,:);
            for j = 1:dim
                distance = abs(toAWinner(1,j) - newAWinner(1,j));
                newAWinner(1,j) = toAWinner(1,j) + (2*rand()-1) * distance;
            end
            newAWFitness=fobj(newAWinner(1,:));
            FEs = FEs + 1;

            %Replace only if improves
            if newAWFitness < AllFitness(i) 
                X(i,:) = newAWinner(1,:);
                AllFitness(i) = newAWFitness;
                aWinners(a,:) = newAWinner(1,:);
            end
        end

        t=t+1;
        Convergence_curve(t)=bestFitness;
    end
    toc
end
