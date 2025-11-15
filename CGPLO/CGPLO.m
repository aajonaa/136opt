function [best_pos,Convergence_curve]=CGPLO(N,MaxFEs,lb,ub,dim,fobj)
% CGPLO: Centroid-Generative Particle Light Optimization
% This algorithm enhances PLO with clustering-based guidance and 
% a population generation mechanism inspired by multi-center concepts.

% --- Internal Algorithm Parameters ---
k_clusters_val = max(2, min(10, ceil(0.1 * N))); 
clustering_freq_val = 10; % Perform clustering every 10 iterations
pg_trigger_probability_val = 0.5; % Probability to trigger PG when conditions met
pg_replacement_fraction_val = 0.1; % Fraction of N particles PG will replace
% --- End Internal Algorithm Parameters ---

%% Initialization
FEs = 0; % Function Evaluations counter
it = 1;  % Iteration counter

AllFitness=inf*ones(N,1); % Stores fitness of all particles
X=initialization(N,dim,ub,lb); % Initialize particle positions

cluster_centroids = []; % Stores centroids of clusters
particle_clusters = []; % Stores cluster assignment for each particle
actual_k = 0; % Actual number of clusters found by k-means

% Evaluate initial population
for i=1:N
    AllFitness(i)=fobj(X(i,:));
    FEs=FEs+1;
end

% Sort initial population
[AllFitness, SortOrder]=sort(AllFitness);
X=X(SortOrder,:);

Bestpos=X(1,:); % Best position found so far
bestFitness=AllFitness(1); % Fitness of the best position

% Pre-allocate convergence curve for efficiency
% The size is an estimate; it will be trimmed at the end.
convergence_curve_max_len = ceil(MaxFEs / N) + 200; % Generous pre-allocation
Convergence_curve=zeros(1, convergence_curve_max_len); 
Convergence_curve(it)=bestFitness;

temp_newX_aurora = zeros(N,dim); % Temporary storage for particles after Aurora walk

%% Main loop
while FEs < MaxFEs % Use < to ensure MaxFEs is not exceeded within the loop's FEs increments
    % --- Innovation 1: Clustering (Periodically) ---
    if mod(it, clustering_freq_val) == 0 || it == 1
        if N >= k_clusters_val % Ensure enough particles to form k_clusters_val
            try
                opts = statset('MaxIter', 50, 'Display', 'off'); % Kmeans options
                % Request k_clusters_val. actual_k will be determined by k-means.
                [~, cluster_centroids_temp, ~, ~, ~] = kmeans(X, k_clusters_val, ...
                    'Distance', 'sqeuclidean', ...
                    'Replicates', 3, ... % Run k-means multiple times to find a good solution
                    'Options', opts, ...
                    'EmptyAction', 'drop'); % Drop clusters if they become empty

                % Filter out any NaN centroids that might result from 'drop' or other issues
                valid_centroids_idx = ~all(isnan(cluster_centroids_temp),2);
                cluster_centroids = cluster_centroids_temp(valid_centroids_idx,:);
                actual_k = size(cluster_centroids, 1);

                if actual_k > 0
                    % Re-assign particle_clusters based on the final valid_centroids
                    % This is important if 'drop' resulted in fewer centroids
                    temp_particle_assignments = zeros(N,1);
                    for p_idx = 1:N
                        % Calculate squared Euclidean distances to all valid centroids
                        dists_to_centroids = sum(bsxfun(@minus, X(p_idx,:), cluster_centroids).^2, 2);
                        [~, temp_particle_assignments(p_idx)] = min(dists_to_centroids);
                    end
                    particle_clusters = temp_particle_assignments;
                else
                    particle_clusters = []; % No valid clusters found
                    % actual_k is already 0
                end
            catch ME
                % disp(['KMeans warning/error: ', ME.message, '. Using global mean for guidance this iteration.']);
                actual_k = 0; 
                cluster_centroids = []; 
                particle_clusters = [];
            end
        else
            % Not enough particles to perform clustering as configured
            % disp(['N < k_clusters_val (', num2str(N), ' < ', num2str(k_clusters_val), '), skipping k-means.']);
            actual_k = 0; 
            cluster_centroids = []; 
            particle_clusters = [];
        end
    end
    % --- End Clustering ---

    X_mean_global=mean(X,1); % Global mean of the population

    % Adaptive weights for PLO
    w1=tansig((FEs/MaxFEs)^4); 
    w2=exp(-(2*FEs/MaxFEs)^3);

    % --- Aurora oval walk (Exploration/Exploitation) ---
    for i=1:N
        if FEs >= MaxFEs; break; end % Check FEs budget before generating new solution

        a=rand()/2+1; % Parameter for LS, range [1, 1.5]
        exp_val = exp((1-a)/100*FEs); % Local search strength, decays as FEs increases
        LS_term_contribution = w1 * exp_val;

        % Determine guidance source for Global Search (GS)
        GS_guidance_source = X_mean_global - X(i,:); % Default: guided by global mean

        % If clustering is active and particle belongs to a valid cluster
        if actual_k > 0 && ~isempty(particle_clusters) && particle_clusters(i) >= 1 && particle_clusters(i) <= actual_k
            X_mean_local = cluster_centroids(particle_clusters(i), :);
            beta_guidance = rand(); % Blending factor
            % Blend local cluster centroid guidance with global best guidance
            GS_guidance_source = beta_guidance * (X_mean_local - X(i,:)) + (1-beta_guidance) * (Bestpos - X(i,:));
        end

        % Levy flight component for GS
        GS_levy_component = Levy(dim).*(GS_guidance_source + (lb+rand(1,dim).*(ub-lb))/2);
        GS_term_contribution = w2 * GS_levy_component; 

        % Calculate step vector and generate new candidate position
        step_vector = (LS_term_contribution * ones(1,dim) + GS_term_contribution);
        newX_i = X(i, :) + step_vector .* rand(1, dim); % Element-wise random scaling of step

        % Boundary checking and correction
        Flag4ub_i=newX_i>ub; 
        Flag4lb_i=newX_i<lb;
        newX_i=(newX_i.*(~(Flag4ub_i+Flag4lb_i)))+ub.*Flag4ub_i+lb.*Flag4lb_i;
        temp_newX_aurora(i,:) = newX_i;
    end
    if FEs >= MaxFEs; break; end % Check FEs after Aurora walk loop

    current_newX = temp_newX_aurora; 

    % --- Particle collision ---
    E =sqrt(FEs/MaxFEs); % Collision probability factor
    A=randperm(N); % Permutation for selecting collision partner
    for i=1:N
        if FEs >= MaxFEs; break; end % Check FEs budget

        collided_X_i = current_newX(i,:); 
        for j_dim=1:dim 
            if (rand<0.05) && (rand<E) % Probability for collision
                 collided_X_i(j_dim) = collided_X_i(j_dim)+sin(rand*pi)*(collided_X_i(j_dim)-current_newX(A(i),j_dim));
            end
        end
        % Boundary checking for collided particle
        Flag4ub_coll=collided_X_i>ub; 
        Flag4lb_coll=collided_X_i<lb;
        collided_X_i=(collided_X_i.*(~(Flag4ub_coll+Flag4lb_coll)))+ub.*Flag4ub_coll+lb.*Flag4lb_coll;
        current_newX(i,:) = collided_X_i; 
    end
    if FEs >= MaxFEs; break; end % Check FEs after collision loop

    % --- Innovation 2: Population Generator (PG) based on cluster centers ---
    % Trigger PG if clustering yielded at least 2 centers and random chance passes
    if actual_k >= 2 && rand() < pg_trigger_probability_val 
        num_pg_particles = ceil(N * pg_replacement_fraction_val); 

        % Replace random particles. Ensure num_pg_particles is not zero.
        if num_pg_particles > 0
            pg_indices_to_replace = randperm(N, num_pg_particles); 
            for k_pg = 1:num_pg_particles
                if FEs >= MaxFEs; break; end % Check FEs budget

                idx_to_replace = pg_indices_to_replace(k_pg);

                % Select two distinct random cluster centroids
                center_indices = randperm(actual_k, 2);
                center1 = cluster_centroids(center_indices(1), :);
                center2 = cluster_centroids(center_indices(2), :);

                dis_vec = center2 - center1; % Vector difference between centers
                rand_scalar = (rand() * 2) - 1; % Random scalar between -1 and 1
                pg_new_particle = center1 + rand_scalar * dis_vec; % Generate new particle

                % Boundary checking for PG particle
                Flag4ub_pg=pg_new_particle>ub; 
                Flag4lb_pg=pg_new_particle<lb;
                pg_new_particle=(pg_new_particle.*(~(Flag4ub_pg+Flag4lb_pg)))+ub.*Flag4ub_pg+lb.*Flag4lb_pg;
                current_newX(idx_to_replace,:) = pg_new_particle;
            end
        end
    end
    if FEs >= MaxFEs; break; end % Check FEs after PG
    % --- End Innovation 2 ---

    % Evaluate new candidate solutions and update population
    newFitness_eval = zeros(N,1); 
    for i=1:N
        if FEs >= MaxFEs; break; end % Check FEs budget before evaluation

        newFitness_eval(i)=fobj(current_newX(i,:)); 
        FEs=FEs+1; % Increment FEs after each evaluation

        if newFitness_eval(i)<AllFitness(i)
            X(i,:)=current_newX(i,:);
            AllFitness(i)=newFitness_eval(i);
        end
    end
    if FEs >= MaxFEs; break; end % Final check for the iteration

    % Sort population based on new fitness values
    [AllFitness, SortOrder]=sort(AllFitness);
    X=X(SortOrder,:); 

    % Update global best
    if AllFitness(1)<bestFitness
        Bestpos=X(1,:);
        bestFitness=AllFitness(1);
    end

    it = it + 1; % Increment iteration counter
    % Record best fitness for convergence curve
    Convergence_curve(it)=bestFitness;
end % End main loop
best_pos=Bestpos;
Convergence_curve = Convergence_curve(1:it);
end % End of CGPLO function


%% Levy flight function (Helper function)
% Generates Levy flight step lengths
% d: dimension of the problem
function o=Levy(d)
beta=1.5; % Levy exponent
% Calculate sigma for Levy distribution
sigma_num = gamma(1+beta)*sin(pi*beta/2);
sigma_den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
sigma=(sigma_num/sigma_den)^(1/beta);

% Generate step vector
u=randn(1,d)*sigma;
v=randn(1,d);
step=u./(abs(v).^(1/beta));
o=0.01*step; % Scaling factor for Levy steps (can be tuned)
% The 0.01 scaling factor is a common heuristic to control step size.
% It might need adjustment based on the problem's landscape.
end