% Nenavath, H., et al. (2018). "A synergy of the sine-cosine algorithm and particle swarm optimizer for improved global optimization and object tracking." Swarm and Evolutionary Computation.

function [SCAgbest, cg_curve]=Hybrid_SCA_PSO(noP,MaxFEs,lb,ub,dim,fobj)
%A synergy of the sine-cosine algorithm and particle swarm optimizer for improved global optimization and object tracking
fes = 0;
%PSO Infotmation
Vmax=4;
wMax=0.9;
wMin=0.4;
c1=2;
c2=2;
PSOMaxIter = 50;

X = initialization(noP,dim,ub,lb);

SCAgbest = zeros(1,dim);
SCAgbestFit = inf;
pbest = zeros(noP,dim);
pbestFit = zeros(noP,1);
for i = 1 : noP
    pbestFit(i) = inf;
end
v = zeros(1,dim);


% main loop
t = 1;
while ( fes <= MaxFEs )
    
    % Eq. (3.4)
    a = 2;
    r1=a-fes*((a)/MaxFEs); % r1 decreases linearly from a to 0
    
    % calculating the fitness of all particles
    for i = 1 : noP

        Fit = fobj(X(i,:));
        fes = fes + 1;

        if (Fit < pbestFit(i))
            pbest(i,:) = X(i,:);
            pbestFit(i) = Fit;

            if (Fit < SCAgbestFit)
                SCAgbest = X(i,:);
                SCAgbestFit = Fit;
            end
        end
    end
    
    for i = 1 : noP
        for k = 1 : dim
            % Update r2, r3, and r4 for Eq. (3.3)
            r2=(2*pi)*rand();
            r3=2*rand;
            r4=rand();

            % Eq. (3.3)
            if r4<0.5
                % Eq. (3.1)
                X(i,:)= X(i,:)+(r1*sin(r2)*abs(r3*SCAgbest(k)-X(i,:)));
            else
                % Eq. (3.2)
                X(i,:)= X(i,:)+(r1*cos(r2)*abs(r3*SCAgbest(k)-X(i,:)));
            end
        end

        % Return back the particles that go beyond the boundaries of the search
        % space
        Flag4ub=X(i,:)>ub;
        Flag4lb=X(i,:)<lb;
        X(i,:)=(X(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
    end
    
    x = pbest;
    PSOgbest = SCAgbest;
    PSOgbestFit = SCAgbestFit;
    
    % PSO loop
    tPSO = 1;
    while ( tPSO <= PSOMaxIter )
        
        %Update the W of PSO
        w=wMax-tPSO*((wMax-wMin)/PSOMaxIter);
        
        % calculating the fitness of particles
        for i = 1 : noP
            Fit = fobj(x(i,:));
            fes = fes + 1;
            
            if (Fit < pbestFit(i))
                pbest(i,:) = x(i,:);
                pbestFit(i) = Fit;

                if (Fit < PSOgbestFit)
                    PSOgbestFit = Fit;
                    PSOgbest = x(i,:);
                end
            end
        end
            
        for i = 1 : noP
            for k = 1 : dim
                v(k)=w*v(k)+c1*rand()*(pbest(i,k)-x(i,k))+c2*rand()*(SCAgbest(k)-x(i,k));

                if(v(k)>Vmax)
                    v(k)=Vmax;
                end
                if(v(k)<-Vmax)
                    v(k)=-Vmax;
                end            
                x(i,k)=x(i,k)+v(k);
            end
        end
        
        % Return back the particles that go beyond the boundaries of the search
        % space
        Flag4ub=x(i,:)>ub;
        Flag4lb=x(i,:)<lb;
        x(i,:)=(x(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        tPSO = tPSO + 1;
    end
    
    cg_curve(t) = PSOgbestFit;
    SCAgbest = PSOgbest;
    SCAgbestFit = PSOgbestFit;
    t = t + 1;
end
end

