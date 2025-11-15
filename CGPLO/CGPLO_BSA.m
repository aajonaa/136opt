function [best_pos,Convergence_curve]=CGPLO_BSA(N,MaxFEs,lb,ub,dim,fobj)
% CGPLO_BSA: Centroid-Generative Particle Light Optimization with BSA-inspired operator
% Enhances CGPLO by adding a recombination mechanism from BSA.

% --- CGPLO Internal Algorithm Parameters ---
k_clusters_val = max(2, min(10, ceil(0.1 * N)));
clustering_freq_val = 10;
pg_trigger_probability_val = 0.5;
pg_replacement_fraction_val = 0.1;
% --- End CGPLO Internal Algorithm Parameters ---

% --- BSA-inspired Operator Internal Parameters ---
DIM_RATE_bsa = 1; % As in the provided BSA example
history_update_prob_bsa = 0.5; % Probability to update historyX_cgplo
map_strategy_prob_bsa = 0.5; % Probability to choose the first map generation strategy
% --- End BSA-inspired Operator Internal Parameters ---

%% Initialization
FEs = 0;
it = 1;
AllFitness=inf*ones(N,1);
X=initialization(N,dim,ub,lb);

% Initialize historical population for BSA-inspired operator
historyX_cgplo = initialization(N,dim,ub,lb); % Or historyX_cgplo = X; after first eval

cluster_centroids = [];
particle_clusters = [];
actual_k = 0;

for i=1:N
    if FEs >= MaxFEs; break; end
    AllFitness(i)=fobj(X(i,:));
    FEs=FEs+1;
end

Bestpos = [];
bestFitness = inf;
if N > 0 && FEs > 0
    [AllFitness, SortOrder]=sort(AllFitness);
    X=X(SortOrder,:);
    Bestpos=X(1,:);
    bestFitness=AllFitness(1);
    historyX_cgplo = X; % Initialize historyX with the sorted initial population
elseif N > 0
    Bestpos = X(1,:);
    historyX_cgplo = X;
end

convergence_curve_max_len = ceil(MaxFEs / N) + 200;
if convergence_curve_max_len < 1; convergence_curve_max_len = 1; end
Convergence_curve=zeros(1, convergence_curve_max_len);
Convergence_curve(it)=bestFitness;

temp_newX_aurora = zeros(N,dim);

%% Main loop
while FEs < MaxFEs
    if N == 0; break; end

    % --- BSA-inspired: SELECTION-I (Update and Shuffle Historical Population) ---
    if rand < history_update_prob_bsa
        historyX_cgplo = X; % Update historical population from current population X
    end
    historyX_cgplo = historyX_cgplo(randperm(N),:); % Shuffle
    % --- End BSA-inspired SELECTION-I ---

    % --- CGPLO: Clustering (Periodically) ---
    if mod(it, clustering_freq_val) == 0 || it == 1
        if N >= k_clusters_val
            try
                opts = statset('MaxIter', 50, 'Display', 'off');
                [~, cluster_centroids_temp] = kmeans(X, k_clusters_val, ...
                    'Distance', 'sqeuclidean', 'Replicates', 3, ...
                    'Options', opts, 'EmptyAction', 'drop');
                valid_centroids_idx = ~all(isnan(cluster_centroids_temp),2);
                cluster_centroids = cluster_centroids_temp(valid_centroids_idx,:);
                actual_k = size(cluster_centroids, 1);
                if actual_k > 0
                    temp_particle_assignments = zeros(N,1);
                    for p_idx = 1:N
                        dists_to_centroids = sum(bsxfun(@minus, X(p_idx,:), cluster_centroids).^2, 2);
                        [~, temp_particle_assignments(p_idx)] = min(dists_to_centroids);
                    end
                    particle_clusters = temp_particle_assignments;
                else particle_clusters = []; end
            catch ME
                actual_k = 0; cluster_centroids = []; particle_clusters = [];
            end
        else actual_k = 0; cluster_centroids = []; particle_clusters = []; end
    end
    % --- End CGPLO Clustering ---

    X_mean_global=mean(X,1);
    w1=tanh((FEs/MaxFEs)^4);
    w2=exp(-(2*FEs/MaxFEs)^3);

    % --- CGPLO: Aurora oval walk ---
    for i=1:N
        if FEs >= MaxFEs; break; end
        a=rand()/2+1; exp_val = exp((1-a)/100*FEs);
        LS_term_contribution = w1 * exp_val;
        GS_guidance_source = X_mean_global - X(i,:);
        if actual_k > 0 && ~isempty(particle_clusters) && particle_clusters(i) >= 1 && particle_clusters(i) <= actual_k
            X_mean_local = cluster_centroids(particle_clusters(i), :);
            beta_guidance = rand();
            GS_guidance_source = beta_guidance * (X_mean_local - X(i,:)) + (1-beta_guidance) * (Bestpos - X(i,:));
        end
        GS_levy_component = Levy(dim).*(GS_guidance_source + (lb+rand(1,dim).*(ub-lb))/2);
        GS_term_contribution = w2 * GS_levy_component;
        step_vector = (LS_term_contribution * ones(1,dim) + GS_term_contribution);
        newX_i = X(i, :) + step_vector .* rand(1, dim);
        Flag4ub_i=newX_i>ub; Flag4lb_i=newX_i<lb;
        newX_i=(newX_i.*(~(Flag4ub_i+Flag4lb_i)))+ub.*Flag4ub_i+lb.*Flag4lb_i;
        temp_newX_aurora(i,:) = newX_i;
    end
    if FEs >= MaxFEs; break; end
    current_newX = temp_newX_aurora;
    % --- End CGPLO Aurora oval walk ---

    % --- CGPLO: Particle collision ---
    E =sqrt(FEs/MaxFEs); A=randperm(N);
    for i=1:N
        if FEs >= MaxFEs; break; end
        collided_X_i = current_newX(i,:);
        for j_dim=1:dim
            if (rand<0.05) && (rand<E)
                 collided_X_i(j_dim) = collided_X_i(j_dim)+sin(rand*pi)*(collided_X_i(j_dim)-current_newX(A(i),j_dim));
            end
        end
        Flag4ub_coll=collided_X_i>ub; Flag4lb_coll=collided_X_i<lb;
        collided_X_i=(collided_X_i.*(~(Flag4ub_coll+Flag4lb_coll)))+ub.*Flag4ub_coll+lb.*Flag4lb_coll;
        current_newX(i,:) = collided_X_i;
    end
    if FEs >= MaxFEs; break; end
    % --- End CGPLO Particle collision ---

    % --- CGPLO: Population Generator (PG) ---
    if actual_k >= 2 && rand() < pg_trigger_probability_val
        num_pg_particles = ceil(N * pg_replacement_fraction_val);
        if num_pg_particles > 0
            actual_num_to_replace = min(N, num_pg_particles);
            if actual_num_to_replace > 0
                pg_indices_to_replace = randperm(N, actual_num_to_replace);
                for k_pg_loop_idx = 1:length(pg_indices_to_replace)
                    if FEs >= MaxFEs; break; end
                    idx_to_replace = pg_indices_to_replace(k_pg_loop_idx);
                    center_indices = randperm(actual_k, 2);
                    center1 = cluster_centroids(center_indices(1), :);
                    center2 = cluster_centroids(center_indices(2), :);
                    dis_vec = center2 - center1; rand_scalar = (rand() * 2) - 1;
                    pg_new_particle = center1 + rand_scalar * dis_vec;
                    Flag4ub_pg=pg_new_particle>ub; Flag4lb_pg=pg_new_particle<lb;
                    pg_new_particle=(pg_new_particle.*(~(Flag4ub_pg+Flag4lb_pg)))+ub.*Flag4ub_pg+lb.*Flag4lb_pg;
                    current_newX(idx_to_replace,:) = pg_new_particle;
                end
            end
        end
    end
    if FEs >= MaxFEs; break; end
    % --- End CGPLO Population Generator ---

    % --- BSA-inspired: RECOMBINATION (MUTATION+CROSSOVER) ---
    F_bsa = get_scale_factor_bsa(); % Get scaling factor
    map_bsa = zeros(N,dim);      % Initialize map
    if rand < map_strategy_prob_bsa
        for i=1:N
            u_bsa=randperm(dim);
            % Ensure ceil argument is not zero if dim is small or rand is very small
            num_dims_to_mutate = ceil(DIM_RATE_bsa*rand*dim);
            if num_dims_to_mutate == 0; num_dims_to_mutate = 1; end % Mutate at least one dim
            map_bsa(i,u_bsa(1:min(dim, num_dims_to_mutate))) = 1;
        end
    else
        for i=1:N
            map_bsa(i,randi(dim)) = 1;
        end
    end
    % Apply BSA recombination to current_newX using historyX_cgplo
    current_newX = current_newX + (map_bsa.*F_bsa).*(historyX_cgplo - current_newX);
    current_newX = BoundaryControl_bsa(current_newX,lb,ub); % Apply BSA boundary control
    % --- End BSA-inspired RECOMBINATION ---

    % --- CGPLO: Evaluate new candidate solutions and update population ---
    newFitness_eval = zeros(N,1);
    for i=1:N
        if FEs >= MaxFEs; break; end
        newFitness_eval(i)=fobj(current_newX(i,:));
        FEs=FEs+1;
        if newFitness_eval(i)<AllFitness(i)
            X(i,:)=current_newX(i,:);
            AllFitness(i)=newFitness_eval(i);
        end
    end
    if FEs >= MaxFEs; break; end

    [AllFitness, SortOrder]=sort(AllFitness);
    X=X(SortOrder,:);
    if AllFitness(1)<bestFitness
        Bestpos=X(1,:);
        bestFitness=AllFitness(1);
    end

    it = it + 1;
    if it <= convergence_curve_max_len
        Convergence_curve(it)=bestFitness;
    else
        Convergence_curve = [Convergence_curve, bestFitness];
        convergence_curve_max_len = length(Convergence_curve);
    end
end % End main loop

best_pos=Bestpos;
if it > 0 && it <= length(Convergence_curve)
    Convergence_curve = Convergence_curve(1:it);
elseif it == 0 && ~isempty(Convergence_curve)
    Convergence_curve = Convergence_curve(1);
else
    Convergence_curve = bestFitness;
end
if isempty(Convergence_curve) && isinf(bestFitness)
    Convergence_curve = inf;
elseif isempty(Convergence_curve)
    Convergence_curve = bestFitness;
end
end % End of CGPLO_BSA function

%% Levy flight function (Helper function for CGPLO)
function o=Levy(d)
beta=1.5; sigma_num = gamma(1+beta)*sin(pi*beta/2);
sigma_den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
sigma=(sigma_num/sigma_den)^(1/beta);
u=randn(1,d)*sigma; v=randn(1,d); step=u./(abs(v).^(1/beta));
o=0.01*step;
end

%% initialization function (Helper function for CGPLO & BSA history)
function X=initialization(N,dim,ub,lb)
    if isscalar(ub); ub_vec = repmat(ub,1,dim); else ub_vec = ub; end
    if isscalar(lb); lb_vec = repmat(lb,1,dim); else lb_vec = lb; end
    if size(lb_vec,1)>1; lb_vec=lb_vec'; end
    if size(ub_vec,1)>1; ub_vec=ub_vec'; end
    X = zeros(N,dim);
    if N > 0
        for i=1:N
            X(i,:) = lb_vec + rand(1,dim).*(ub_vec-lb_vec);
        end
    end
end

%% BSA-inspired Helper Functions
% F factor for BSA-inspired operator
function F=get_scale_factor_bsa()
    F=3*randn; % STANDARD brownian-walk (as per BSA example)
    % Other options from BSA example could be used here if desired
    % F=4*randg;
    % F=lognrnd(rand,5*rand);
    % F=1/normrnd(0,5);
    % F=1./gamrnd(1,0.5);
end

% Boundary control for BSA-inspired operator
function X=BoundaryControl_bsa(X,lb,ub)
[N,dim]=size(X);
% Ensure lb and ub are row vectors of length dim for comparison
if isscalar(lb); lb_dim = repmat(lb,1,dim); else lb_dim = lb; end
if isscalar(ub); ub_dim = repmat(ub,1,dim); else ub_dim = ub; end
if size(lb_dim,1)>1; lb_dim=lb_dim'; end
if size(ub_dim,1)>1; ub_dim=ub_dim'; end

for i=1:N
    for j=1:dim
        k_boundary_strat = rand < 0.5; % Corresponds to rand < rand in BSA example
        if X(i,j)<lb_dim(j)
            if k_boundary_strat, X(i,j)=lb_dim(j);
            else X(i,j)=rand*(ub_dim(j)-lb_dim(j))+lb_dim(j);
            end
        end
        if X(i,j)>ub_dim(j)
            if k_boundary_strat, X(i,j)=ub_dim(j);
            else X(i,j)=rand*(ub_dim(j)-lb_dim(j))+lb_dim(j);
            end
        end
    end
end
end