% Issa, M., et al. (2018). "ASCA-PSO: Adaptive sine cosine optimization algorithm integrated with particle swarm for pairwise local sequence alignment." Expert Systems with Applications 99: 56-70.

function [ygbest, cg_curve]=ASCA_PSO(noP,MaxFEs,lb,ub,dim,fobj)
if (size(ub, 2)~=1)
    dim = size(ub, 2);
else
    lb = lb*ones(1,dim);
    ub = ub*ones(1,dim);
end
%Adaptive sine cosine integrated with particle swarm (ASCA-PSO)
fes = 0;
M = 4;
N = 9;
if ( noP < M*N )
    noP = M*N;
end
% population size = M * N + M

%PSO Infotmation
Vmax=6;
wMax=0.9;
wMin=0.2;
c1=2;
c2=2;

X = initialization(noP,dim,ub,lb);
x = zeros(M,N,dim);
xFit = zeros(M,N);
y = zeros(M,dim);
yFit = zeros(M,1);
for i = 1 : M
    yFit(i) = inf;
end
groupBestFit = zeros(M,1);
for i = 1 : M
    groupBestFit(i) = inf;
end
ygbest = zeros(1,dim);
ygbestFit = inf;
ypbest = zeros(M,dim);
ypbestFit = zeros(M,1);
for i = 1 : M
    ypbestFit(i) = inf;
end
v = zeros(1,dim);

%initialize x
k = 1;
for i = 1 : M
    for j = 1 : N
        x(i,j,:) = X(k,:);
        k = k + 1;
    end
end

%initialize y
k = 1;
for i = 1 : M
    y(i,:) = X(k,:);
    k = k + 1;
end

% First time calculating the fitness of all particles
for i = 1 : M
    for j = 1 : N
        tmpX = zeros(1, dim);
        for p = 1 : dim
            tmpX(p) = x(i, j, p);
        end
        xFit(i,j) = fobj(tmpX);
        fes = fes + 1;
        
        if (xFit(i,j) < groupBestFit(i))
            groupBestFit(i) = xFit(i,j);

            if (groupBestFit(i) < yFit(i))
                y(i,:) = x(i,j,:);
                yFit(i) = xFit(i,j);

                if (yFit(i) < ygbestFit)
                    ygbest = y(i,:);
                    ygbestFit = yFit(i);
                end
            end
        end
    end
end

for i = 1 : M
    yFit(i) = fobj(y(i,:));
    fes = fes + 1;
    
    if (yFit(i) < ypbestFit(i))
        ypbest(i,:) = y(i,:);
        ypbestFit(i) = yFit(i);
        
        if (yFit(i) < ygbestFit)
            ygbest = y(i,:);
            ygbestFit = yFit(i);
        end
    end
end


% main loop
t = 1;
while ( fes <= MaxFEs )
    
    % Eq. (3.4)
    a = 2;
    r1=a-fes*((a)/MaxFEs); % r1 decreases linearly from a to 0
    
    %Update the W of PSO
    w=wMax-fes*((wMax-wMin)/MaxFEs);
    
    for i = 1 : M
        for j = 1 : N
            for k = 1 : dim
                % Update r2, r3, and r4 for Eq. (3.3)
                r2=(2*pi)*rand();
                r3=2*rand;
                r4=rand();
                
                % Eq. (3.3)
                if r4<0.5
                    % Eq. (3.1)
                    x(i,j,k)= x(i,j,k)+(r1*sin(r2)*abs(r3*y(i,k)-x(i,j,k)));
                else
                    % Eq. (3.2)
                    x(i,j,k)= x(i,j,k)+(r1*cos(r2)*abs(r3*y(i,k)-x(i,j,k)));
                end
            end
            
            % Return back the particles that go beyond the boundaries of the search
            % space
            x1 = zeros(1, dim);
            for k = 1 : dim
                Flag4ub(k)=x(i, j, k)>ub(k);
                Flag4lb(k)=x(i, j, k)<lb(k);
                x1(k) = x(i, j, k);
            end
            x1=(x1.*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
            x(i, j, :) = x1;
            
            % calculating the fitness of particles
            tmpX = zeros(1, dim);
            for p = 1 : dim
                tmpX(p) = x(i, j, p);
            end
            xFit(i,j) = fobj(tmpX);
            fes = fes + 1;

            if (xFit(i,j) < groupBestFit(i))
                groupBestFit(i) = xFit(i,j);
                
                if (groupBestFit(i) < yFit(i))
                    y(i,:) = x(i,j,:);
                    yFit(i) = xFit(i,j);

                    if (yFit(i) < ygbestFit)
                        ygbestFit = yFit(i);
                        ygbest = y(i,:);
                    end
                end
            end
        end
        
        % PSO
        for k = 1 : dim
            v(k)=w*v(k)+c1*rand()*(ypbest(i,k)-y(i,k))+c2*rand()*(ygbest(k)-y(i,k));

            if(v(k)>Vmax)
                v(k)=Vmax;
            end
            if(v(k)<-Vmax)
                v(k)=-Vmax;
            end            
            y(i,k)=y(i,k)+v(k);
        end
        
        % Return back the particles that go beyond the boundaries of the search
        % space
        Flag4ub=y(i,:)>ub;
        Flag4lb=y(i,:)<lb;
        y(i,:)=(y(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
        
        % calculating the fitness of particles
        yFit(i) = fobj(y(i,:));
        fes = fes + 1;
    
        if (yFit(i) < ypbestFit(i))
            ypbest(i,:) = y(i,:);
            ypbestFit(i) = yFit(i);

            if (yFit(i) < ygbestFit)
                ygbest = y(i,:);
                ygbestFit = yFit(i);
            end
        end
    end
    
    cg_curve(t) = ygbestFit;
    t = t + 1;
end
end

