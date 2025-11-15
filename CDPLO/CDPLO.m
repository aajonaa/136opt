% function [best_pos,Convergence_curve]=CDPLO(N,MaxFEs,lb,ub,dim,fobj)
% % CGPLO_DR_v4: PLO with Modular Innovations
% % Innovation 1: Additive Cluster-Guided Aurora Oval Walk (ACG-AW)
% % Innovation 2: BSA-Inspired Differential Recombination (BSA-DR)
% 
% % --- Control Flags for Innovations ---
% use_innovation1 = true; % Set to false to disable Additive Cluster-Guided Aurora Walk
% use_innovation2 = true; % Set to false to disable BSA-Inspired Differential Recombination
% % --- End Control Flags ---
% 
% % --- ACG-AW Internal Parameters (Used if use_innovation1 is true) ---
% k_clusters_val = max(2, min(10, ceil(0.1 * N)));
% clustering_freq_val = 10;
% cluster_influence_factor = 0.2; % Strength of the additive cluster influence
% % --- End ACG-AW Internal Parameters ---
% 
% % --- BSA-DR Internal Parameters (Used if use_innovation2 is true) ---
% history_update_prob = 0.5; 
% DIM_RATE_bsa_dr = 1;       
% map_strategy_prob_bsa_dr = 0.5; 
% % --- End BSA-DR Internal Parameters ---
% 
% %% Initialization
% FEs = 0;
% it = 1;
% AllFitness=inf*ones(N,1);
% X=initialization(N,dim,ub,lb);
% 
% historyX = []; % Initialize empty, will be populated if use_innovation2
% if use_innovation2
%     historyX = initialization(N,dim,ub,lb); 
% end
% 
% cluster_centroids = [];
% particle_clusters = [];
% actual_k = 0;
% 
% for i=1:N
%     if FEs >= MaxFEs; break; end
%     AllFitness(i)=fobj(X(i,:));
%     FEs=FEs+1;
% end
% 
% Bestpos = [];
% bestFitness = inf;
% if N > 0 && FEs > 0
%     [AllFitness, SortOrder]=sort(AllFitness);
%     %% -- Algorithm Operator -- (Population sorting is an update to X's order)
%     X=X(SortOrder,:);
%     Bestpos=X(1,:);
%     bestFitness=AllFitness(1);
%     if use_innovation2
%         %% -- Algorithm Operator -- (History population update)
%         historyX = X; 
%     end
% elseif N > 0
%     Bestpos = X(1,:); 
%     if use_innovation2
%         %% -- Algorithm Operator -- (History population update)
%         historyX = X;     
%     end
% end
% 
% convergence_curve_max_len = ceil(MaxFEs / N) + 200;
% if convergence_curve_max_len < 1; convergence_curve_max_len = 1; end
% Convergence_curve=zeros(1, convergence_curve_max_len);
% Convergence_curve(it)=bestFitness;
% 
% temp_newX_aurora = zeros(N,dim);
% 
% %% Main loop
% while FEs < MaxFEs
%     if N == 0; break; end
% 
%     if use_innovation2
%         % --- BSA-DR: Update and Shuffle Historical Population ---
%         if rand < history_update_prob
%             %% -- Algorithm Operator -- (History population update)
%             historyX = X; 
%         end
%         %% -- Algorithm Operator -- (History population shuffling)
%         historyX = historyX(randperm(N),:); 
%         % --- End BSA-DR Historical Population Update ---
%     end
% 
%     if use_innovation1
%         % --- Innovation 1: ACG-AW - Clustering (Periodically) ---
%         if mod(it, clustering_freq_val) == 0 || it == 1
%             if N >= k_clusters_val
%                 try
%                     opts = statset('MaxIter', 50, 'Display', 'off');
%                     [~, cluster_centroids_temp] = kmeans(X, k_clusters_val, ...
%                         'Distance', 'sqeuclidean', 'Replicates', 3, ...
%                         'Options', opts, 'EmptyAction', 'drop');
%                     valid_centroids_idx = ~all(isnan(cluster_centroids_temp),2);
%                     cluster_centroids = cluster_centroids_temp(valid_centroids_idx,:);
%                     actual_k = size(cluster_centroids, 1);
%                     if actual_k > 0
%                         temp_particle_assignments = zeros(N,1);
%                         for p_idx = 1:N
%                             dists_to_centroids = sum(bsxfun(@minus, X(p_idx,:), cluster_centroids).^2, 2);
%                             [~, temp_particle_assignments(p_idx)] = min(dists_to_centroids);
%                         end
%                         particle_clusters = temp_particle_assignments;
%                     else particle_clusters = []; end
%                 catch ME
%                     actual_k = 0; cluster_centroids = []; particle_clusters = [];
%                 end
%             else actual_k = 0; cluster_centroids = []; particle_clusters = []; end
%         end
%         % --- End ACG-AW Clustering ---
%     else
%         actual_k = 0; % Ensure actual_k is 0 if innovation 1 is not used
%         cluster_centroids = [];
%         particle_clusters = [];
%     end
% 
%     X_mean_global=mean(X,1);
%     w1=tansig((FEs/MaxFEs)^4); % Using tansig as in original PLO example
%     w2=exp(-(2*FEs/MaxFEs)^3); 
% 
%     % --- Aurora Oval Walk ---
%     % The modification for Innovation 1 is handled inside this loop
%     for i=1:N
%         if FEs >= MaxFEs; break; end
%         a=rand()/2+1; exp_val = exp((1-a)/100*FEs);
%         LS_term_contribution = w1 * exp_val; % Local search component
% 
%         GS_base_diff_vector = X_mean_global - X(i,:);
%         additive_cluster_influence = zeros(1,dim); % Default to no influence
% 
%         if use_innovation1
%             % --- Innovation 1: ACG-AW - Calculate Additive Cluster Influence ---
%             if actual_k > 0 && ~isempty(particle_clusters) && particle_clusters(i) >= 1 && particle_clusters(i) <= actual_k
%                 additive_cluster_influence = cluster_influence_factor * (cluster_centroids(particle_clusters(i), :) - X(i,:));
%             end
%             % --- End Innovation 1: ACG-AW - Additive Influence Calculation ---
%         end
% 
%         GS_final_guidance_diff_vector = GS_base_diff_vector + additive_cluster_influence;
% 
%         GS_levy_component = Levy(dim).*(GS_final_guidance_diff_vector + (lb+rand(1,dim).*(ub-lb))/2);
%         GS_term_contribution = w2 * GS_levy_component; 
% 
%         step_vector = (LS_term_contribution * ones(1,dim) + GS_term_contribution);
%         %% -- Algorithm Operator --
%         newX_i = X(i, :) + step_vector .* rand(1, dim); 
% 
%         Flag4ub_i=newX_i>ub; Flag4lb_i=newX_i<lb; 
%         %% -- Algorithm Operator -- (Boundary control is a form of update)
%         newX_i=(newX_i.*(~(Flag4ub_i+Flag4lb_i)))+ub.*Flag4ub_i+lb.*Flag4lb_i;
%         temp_newX_aurora(i,:) = newX_i;
%     end
%     if FEs >= MaxFEs; break; end
%     current_newX = temp_newX_aurora; 
%     % --- End Aurora Oval Walk ---
% 
%     % --- Base PLO: Particle Collision ---
%     E =sqrt(FEs/MaxFEs); A=randperm(N);
%     for i=1:N
%         if FEs >= MaxFEs; break; end
%         collided_X_i = current_newX(i,:);
%         for j_dim=1:dim
%             if (rand<0.05) && (rand<E) 
%                  %% -- Algorithm Operator --
%                  collided_X_i(j_dim) = collided_X_i(j_dim)+sin(rand*pi)*(collided_X_i(j_dim)-current_newX(A(i),j_dim));
%             end
%         end
%         Flag4ub_coll=collided_X_i>ub; Flag4lb_coll=collided_X_i<lb; 
%         %% -- Algorithm Operator -- (Boundary control is a form of update)
%         collided_X_i=(collided_X_i.*(~(Flag4ub_coll+Flag4lb_coll)))+ub.*Flag4ub_coll+lb.*Flag4lb_coll;
%         current_newX(i,:) = collided_X_i; 
%     end
%     if FEs >= MaxFEs; break; end
%     % --- End Base PLO Particle Collision ---
% 
%     if use_innovation2
%         % --- Innovation 2: BSA-Inspired Differential Recombination (BSA-DR) ---
%         F_bsa_dr = get_scale_factor_bsa_dr(); 
%         map_bsa_dr = zeros(N,dim);      
% 
%         if rand < map_strategy_prob_bsa_dr
%             for i_map=1:N % Changed loop variable to avoid conflict
%                 u_bsa_dr=randperm(dim);
%                 num_dims_to_mutate = ceil(DIM_RATE_bsa_dr*rand*dim);
%                 if num_dims_to_mutate == 0; num_dims_to_mutate = 1; end
%                 map_bsa_dr(i_map,u_bsa_dr(1:min(dim, num_dims_to_mutate))) = 1;
%             end
%         else
%             for i_map=1:N % Changed loop variable to avoid conflict
%                 map_bsa_dr(i_map,randi(dim)) = 1; 
%             end
%         end
% 
%         differential_term = historyX - current_newX; 
%         %% -- Algorithm Operator --
%         current_newX = current_newX + F_bsa_dr * (map_bsa_dr .* differential_term);
%         %% -- Algorithm Operator -- (Boundary control is a form of update)
%         current_newX = BoundaryControl_bsa_dr(current_newX,lb,ub);
%         % --- End Innovation 2: BSA-DR ---
%     end
% 
%     % --- Evaluation and Selection ---
%     newFitness_eval = zeros(N,1);
%     for i=1:N
%         if FEs >= MaxFEs; break; end
%         newFitness_eval(i)=fobj(current_newX(i,:));
%         FEs=FEs+1;
%         if newFitness_eval(i)<AllFitness(i)
%             %% -- Algorithm Operator -- (Greedy selection updates X)
%             X(i,:)=current_newX(i,:);
%             AllFitness(i)=newFitness_eval(i);
%         end
%     end
%     if FEs >= MaxFEs; break; end
% 
%     [AllFitness, SortOrder]=sort(AllFitness);
%     %% -- Algorithm Operator -- (Population sorting is an update to X's order and content based on fitness)
%     X=X(SortOrder,:); 
% 
%     if N > 0 && FEs > 0 && ~isempty(AllFitness) && AllFitness(1)<bestFitness 
%         Bestpos=X(1,:);
%         bestFitness=AllFitness(1);
%     end
% 
%     it = it + 1;
%     if it <= convergence_curve_max_len
%         Convergence_curve(it)=bestFitness;
%     else
%         Convergence_curve = [Convergence_curve, bestFitness]; 
%         convergence_curve_max_len = length(Convergence_curve);
%     end
% end % End main loop
% 
% best_pos=Bestpos;
% if it > 0 && it <= length(Convergence_curve)
%     Convergence_curve = Convergence_curve(1:it);
% elseif it == 0 && ~isempty(Convergence_curve) 
%     Convergence_curve = Convergence_curve(1);
% else 
%     Convergence_curve = bestFitness;
% end
% if isempty(Convergence_curve) && isinf(bestFitness) 
%     Convergence_curve = inf;
% elseif isempty(Convergence_curve)
%     Convergence_curve = bestFitness;
% end
% end % End of CGPLO_DR_v4 function
% 
% %% Levy flight function (Helper for ACG-AW)
% function o=Levy(d)
% beta=1.5; sigma_num = gamma(1+beta)*sin(pi*beta/2);
% sigma_den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
% sigma=(sigma_num/sigma_den)^(1/beta);
% u=randn(1,d)*sigma; v=randn(1,d); step=u./(abs(v).^(1/beta));
% o=step; % Removed 0.01 scaling to match original PLO behavior
% end
% 
% %% initialization function (Helper for population and historyX)
% function X=initialization(N,dim,ub,lb)
%     if isscalar(ub); ub_vec = repmat(ub,1,dim); else ub_vec = ub; end
%     if isscalar(lb); lb_vec = repmat(lb,1,dim); else lb_vec = lb; end
%     if size(lb_vec,1)>1; lb_vec=lb_vec'; end 
%     if size(ub_vec,1)>1; ub_vec=ub_vec'; end 
%     X = zeros(N,dim);
%     if N > 0
%         for i_init=1:N % Changed loop variable to avoid conflict
%             %% -- Algorithm Operator -- (Initial population generation)
%             X(i_init,:) = lb_vec + rand(1,dim).*(ub_vec-lb_vec);
%         end
%     end
% end
% 
% %% BSA-DR Inspired Helper Functions
% function F=get_scale_factor_bsa_dr()
%     F=3*randn; 
% end
% 
% function X_out=BoundaryControl_bsa_dr(X_in,lb,ub) 
% [N_bc,dim_bc]=size(X_in); % Use different var names to avoid conflict with N, dim
% X_out = X_in; 
% if isscalar(lb); lb_dim = repmat(lb,1,dim_bc); else lb_dim = lb; end
% if isscalar(ub); ub_dim = repmat(ub,1,dim_bc); else ub_dim = ub; end
% if size(lb_dim,1)>1; lb_dim=lb_dim'; end
% if size(ub_dim,1)>1; ub_dim=ub_dim'; end
% 
% for i_bc=1:N_bc % Changed loop variable
%     for j_bc=1:dim_bc % Changed loop variable
%         k_boundary_strat = rand < 0.5; 
%         if X_out(i_bc,j_bc)<lb_dim(j_bc)
%             if k_boundary_strat
%                 %% -- Algorithm Operator -- (Boundary control update)
%                 X_out(i_bc,j_bc)=lb_dim(j_bc); 
%             else
%                 %% -- Algorithm Operator -- (Boundary control update)
%                 X_out(i_bc,j_bc)=rand*(ub_dim(j_bc)-lb_dim(j_bc))+lb_dim(j_bc); 
%             end
%         end
%         if X_out(i_bc,j_bc)>ub_dim(j_bc)
%             if k_boundary_strat
%                 %% -- Algorithm Operator -- (Boundary control update)
%                 X_out(i_bc,j_bc)=ub_dim(j_bc); 
%             else
%                 %% -- Algorithm Operator -- (Boundary control update)
%                 X_out(i_bc,j_bc)=rand*(ub_dim(j_bc)-lb_dim(j_bc))+lb_dim(j_bc); 
%             end
%         end
%     end
% end
% end

function [best_pos,Convergence_curve]=CDPLO(N,MaxFEs,lb,ub,dim,fobj)
% CDPLO: PLO with Additive Cluster-Guided Aurora Oval Walk 
%        and BSA-Inspired Differential Recombination

k_clusters_val = max(2, min(10, ceil(0.1 * N)));
clustering_freq_val = 10;
cluster_influence_factor = 0.2; 

history_update_prob = 0.5; 
DIM_RATE_bsa_dr = 1;       
map_strategy_prob_bsa_dr = 0.5; 

%% Initialization
FEs = 0;
it = 1;
AllFitness=inf*ones(N,1);
X=initialization(N,dim,ub,lb);
historyX = initialization(N,dim,ub,lb); 

cluster_centroids = [];
particle_clusters = [];
actual_k = 0;

for i=1:N
    if FEs >= MaxFEs; break; end
    AllFitness(i)=fobj(X(i,:));
    FEs=FEs+1;
end

[AllFitness, SortOrder]=sort(AllFitness);
X=X(SortOrder,:);
Bestpos=X(1,:);
bestFitness=AllFitness(1);
historyX = X; 


Convergence_curve=[];
Convergence_curve(it)=bestFitness;

temp_newX_aurora = zeros(N,dim);

%% Main loop
while FEs < MaxFEs

    % Update and Shuffle Historical Population
    if rand < history_update_prob
        historyX = X; 
    end
    historyX = historyX(randperm(N),:); 

    % --- Innovation 1: Clustering (Periodically) ---
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
                else 
                    particle_clusters = []; 
                    actual_k = 0; % Ensure actual_k is 0 if no valid centroids
                end
            catch ME
                actual_k = 0; cluster_centroids = []; particle_clusters = [];
            end
        else 
            actual_k = 0; cluster_centroids = []; particle_clusters = []; 
        end
    end
    % --- End Clustering ---

    X_mean_global=mean(X,1);
    w1=tansig((FEs/MaxFEs)^4);
    w2=exp(-(2*FEs/MaxFEs)^3); 

    % --- Innovation 1: Aurora Oval Walk (Additive Cluster Influence) ---
    for i=1:N
        if FEs >= MaxFEs; break; end
        a=rand()/2+1; exp_val = exp((1-a)/100*FEs);
        LS_term_contribution = w1 * exp_val; 
        
        GS_base_diff_vector = X_mean_global - X(i,:);
        additive_cluster_influence = zeros(1,dim); 
        
        % Additive Cluster Influence (Part of Innovation 1)
        if actual_k > 0 && ~isempty(particle_clusters) && particle_clusters(i) >= 1 && particle_clusters(i) <= actual_k
            additive_cluster_influence = cluster_influence_factor * (cluster_centroids(particle_clusters(i), :) - X(i,:));
        end
        
        GS_final_guidance_diff_vector = GS_base_diff_vector + additive_cluster_influence;
        
        GS_levy_component = Levy(dim).*(GS_final_guidance_diff_vector + (lb+rand(1,dim).*(ub-lb))/2);
        GS_term_contribution = w2 * GS_levy_component; 
        
        step_vector = (LS_term_contribution * ones(1,dim) + GS_term_contribution);
        newX_i = X(i, :) + step_vector .* rand(1, dim); 
        
        Flag4ub_i=newX_i>ub; Flag4lb_i=newX_i<lb; 
        newX_i=(newX_i.*(~(Flag4ub_i+Flag4lb_i)))+ub.*Flag4ub_i+lb.*Flag4lb_i;
        temp_newX_aurora(i,:) = newX_i;
    end
    if FEs >= MaxFEs; break; end
    current_newX = temp_newX_aurora; 
    % --- End Innovation 1 ---

    % --- Base PLO: Particle Collision ---
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
    % --- End Base PLO Particle Collision ---

    % --- Innovation 2: Differential Recombination ---
    F_bsa_dr = get_scale_factor_bsa_dr(); 
    map_bsa_dr = zeros(N,dim);      
    
    if rand < map_strategy_prob_bsa_dr
        for i_map=1:N 
            u_bsa_dr=randperm(dim);
            num_dims_to_mutate = ceil(DIM_RATE_bsa_dr*rand*dim);
            if num_dims_to_mutate == 0; num_dims_to_mutate = 1; end
            map_bsa_dr(i_map,u_bsa_dr(1:min(dim, num_dims_to_mutate))) = 1;
        end
    else
        for i_map=1:N 
            map_bsa_dr(i_map,randi(dim)) = 1; 
        end
    end
    
    differential_term = historyX - current_newX; 
    current_newX = current_newX + F_bsa_dr * (map_bsa_dr .* differential_term);
    current_newX = BoundaryControl(current_newX,lb,ub);
    % --- End Innovation 2 ---

    % --- Evaluation and Selection ---
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
    Bestpos=X(1,:);
    bestFitness=AllFitness(1);

    it = it + 1;
    Convergence_curve(it)=bestFitness;
end % End main loop
best_pos=Bestpos;
Convergence_curve = Convergence_curve(1:it);
end % End of CDPLO function

%% Levy flight function (Helper for ACG-AW)
function o=Levy(d)
beta=1.5; sigma_num = gamma(1+beta)*sin(pi*beta/2);
sigma_den = gamma((1+beta)/2)*beta*2^((beta-1)/2);
sigma=(sigma_num/sigma_den)^(1/beta);
u=randn(1,d)*sigma; v=randn(1,d); step=u./(abs(v).^(1/beta));
o=step; % Original PLO Levy step (no 0.01 scaling)
end

function F=get_scale_factor_bsa_dr()
    F=3*randn; 
end

function X=BoundaryControl(X,lb,ub)
    if numel(lb)
        lb = repmat(lb, 1, 30);
        ub = repmat(ub, 1, 30);
    end
    [N,dim]=size(X);
    for i=1:N
        for j=1:dim                
            k=rand<rand; 
            if X(i,j)<lb(j)
                if k, X(i,j)=lb(j); 
                else X(i,j)=rand*(ub(j)-lb(j))+lb(j); 
                end 
            end        
            if X(i,j)>ub(j)
                if k, X(i,j)=ub(j);  
                else
                    X(i,j)=rand*(ub(j)-lb(j))+lb(j); 
                end 
            end
        end
    end
end
