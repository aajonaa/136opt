# import torch
# import torch.nn as nn
# import torch.nn.functional as F
# import torch.optim as optim
# import torch.multiprocessing as mp
# import numpy as np
# import os
# import time
# from collections import deque
# 
# # Set up multiprocessing support
# try:
#     mp.set_start_method('spawn', force=True)
# except RuntimeError:
#     pass  # Already set
# 
# # Define Population & Attention Parameters
# POP_SIZE = 30
# IND_SIZE = 30
# GAMMA = 0.99
# LR = 0.0003
# NUM_PROCESSES = os.cpu_count() - 1 or 1
# MAX_EPISODE_LENGTH = 100
# ENTROPY_BETA = 0.01
# 
# # Shared counter for tracking global steps
# class Counter:
#     def __init__(self):
#         self.val = mp.Value('i', 0)
# 
#     def increment(self):
#         with self.val.get_lock():
#             self.val.value += 1
# 
#     def get(self):
#         with self.val.get_lock():
#             return self.val.value
# 
# global_step_counter = Counter()
# episode_rewards = mp.Queue()
# 
# # Actor-Critic Network with Attention and Layer Normalization
# class ActorCriticWithAttention(nn.Module):
#     def __init__(self, state_dim, action_dim, hidden_dim=128):
#         super(ActorCriticWithAttention, self).__init__()
# 
#         # Shared feature extraction layers
#         self.shared_fc1 = nn.Linear(state_dim, hidden_dim)
#         self.shared_fc2 = nn.Linear(hidden_dim, hidden_dim)
# 
#         # Layer normalization for stability
#         self.layer_norm1 = nn.LayerNorm(hidden_dim)
#         self.layer_norm2 = nn.LayerNorm(hidden_dim)
# 
#         # Self-attention mechanism
#         self.query = nn.Linear(hidden_dim, hidden_dim)
#         self.key = nn.Linear(hidden_dim, hidden_dim)
#         self.value = nn.Linear(hidden_dim, hidden_dim)
# 
#         # Layer normalization for attention
#         self.attention_norm = nn.LayerNorm(hidden_dim)
# 
#         # Actor (policy) head
#         self.actor = nn.Linear(hidden_dim, action_dim)
# 
#         # Critic (value) head
#         self.critic = nn.Linear(hidden_dim, 1)
# 
#     def forward(self, x):
#         # Shared feature extraction with layer normalization
#         x = F.relu(self.layer_norm1(self.shared_fc1(x)))
#         x = F.relu(self.layer_norm2(self.shared_fc2(x)))
# 
#         # Self-attention
#         batch_size = x.size(0)
#         query = self.query(x)
#         key = self.key(x)
#         value = self.value(x)
# 
#         # Compute attention scores with temperature scaling for stability
#         temperature = np.sqrt(key.size(1))
#         if batch_size > 1:
#             attention_scores = torch.bmm(query.unsqueeze(1), key.unsqueeze(2)).squeeze(1) / temperature
#             attention_probs = F.softmax(attention_scores, dim=1)
#             context = attention_probs.unsqueeze(1) * value
#             context = context.sum(dim=1)
#         else:
#             attention_scores = torch.mm(query, key.t()) / temperature
#             attention_probs = F.softmax(attention_scores, dim=1)
#             context = torch.mm(attention_probs, value)
# 
#         # Combine with residual connection and normalization
#         x = self.attention_norm(x + context)
# 
#         # Actor output (action probabilities)
#         logits = self.actor(x)
#         action_probs = F.softmax(logits, dim=-1)
# 
#         # Critic output (state value)
#         state_value = self.critic(x)
# 
#         return action_probs, state_value, attention_probs
# 
# # Worker process for A3C
# class Worker(mp.Process):
#     def __init__(self, global_model, optimizer, global_step_counter, episode_rewards, 
#                  process_idx, fitness_queue, result_queue):
#         super(Worker, self).__init__()
#         self.process_idx = process_idx
#         self.global_model = global_model
#         self.optimizer = optimizer
#         self.global_step_counter = global_step_counter
#         self.episode_rewards = episode_rewards
#         self.fitness_queue = fitness_queue  # Queue to receive external fitness values
#         self.result_queue = result_queue    # Queue to send back modified populations
# 
#         # Local model for this worker
#         self.local_model = ActorCriticWithAttention(IND_SIZE, POP_SIZE)
#         self.local_model.load_state_dict(self.global_model.state_dict())
# 
#         # Episode information
#         self.done = True
#         self.episode_reward = 0
#         self.population = None
#         self.state = None
#         self.fitness_scores = None
# 
#     def run(self):
#         while self.global_step_counter.get() < 10000:  # Limit total steps
#             # Initialize or reset environment
#             if self.done:
#                 # Generate new random population as initial state
#                 self.population = np.random.uniform(-100, 100, (POP_SIZE, IND_SIZE))
# 
#                 # Instead of calculating fitness internally, wait for external fitness
#                 try:
#                     # Non-blocking check for fitness values (timeout after 0.1 seconds)
#                     self.fitness_scores = self.fitness_queue.get(timeout=0.1)
#                 except:
#                     # If no external fitness provided, use simple proxy fitness
#                     self.fitness_scores = -np.sum(np.abs(self.population), axis=1)
# 
#                 # Enhanced state representation: mean, variance, and best individual
#                 pop_mean = self.population.mean(axis=0)
#                 pop_var = self.population.var(axis=0)
#                 best_ind = self.population[np.argmin(self.fitness_scores)]
# 
#                 # Concatenate statistics for richer state representation
#                 state_components = [pop_mean, pop_var, best_ind]
#                 state_vector = np.concatenate(state_components)
# 
#                 self.state = torch.tensor(state_vector, dtype=torch.float32).unsqueeze(0)
#                 self.done = False
#                 self.episode_reward = 0
# 
#             self.local_model.load_state_dict(self.global_model.state_dict())
# 
#             # Storage for trajectory
#             log_probs = []
#             values = []
#             rewards = []
#             entropies = []
# 
#             # Collect trajectory
#             for _ in range(MAX_EPISODE_LENGTH):
#                 # Get action probabilities and state value
#                 action_probs, value, attention_weights = self.local_model(self.state)
# 
#                 # Sample action
#                 dist = torch.distributions.Categorical(action_probs)
#                 action = dist.sample()
# 
#                 # Calculate entropy for exploration encouragement
#                 entropy = -(action_probs * action_probs.log()).sum()
# 
#                 # Convert attention weights to numpy and ensure they're valid
#                 attention_weights = attention_weights.detach().numpy().flatten()
#                 if np.isnan(attention_weights).any() or np.isinf(attention_weights).any():
#                     attention_weights = np.ones(POP_SIZE) / POP_SIZE
# 
#                 # Generate new population through crossover and mutation
#                 new_population = np.zeros_like(self.population)
#                 for i in range(POP_SIZE):
#                     # Select parents using attention weights
#                     parent_indices = np.random.choice(POP_SIZE, 2, p=attention_weights/attention_weights.sum())
#                     parent1, parent2 = self.population[parent_indices[0]], self.population[parent_indices[1]]
# 
#                     # Crossover
#                     crossover_point = np.random.randint(1, IND_SIZE)
#                     child = np.concatenate([parent1[:crossover_point], parent2[crossover_point:]])
# 
#                     # Adaptive mutation based on training progress
#                     mutation_rate = max(0.05, 0.2 * (1 - self.global_step_counter.get() / 10000))
#                     mutation_mask = np.random.random(IND_SIZE) < mutation_rate
#                     mutation_strength = 0.1 * (1 - self.global_step_counter.get() / 10000)
#                     child[mutation_mask] += np.random.normal(0, mutation_strength, size=mutation_mask.sum())
# 
#                     new_population[i] = child
# 
#                 # Try to get external fitness, fall back to proxy if not available
#                 try:
#                     # Send population for evaluation and get fitness back
#                     self.result_queue.put(new_population)
#                     new_fitness_scores = self.fitness_queue.get(timeout=0.1)
#                 except:
#                     # Use proxy fitness if external evaluation is not available
#                     new_fitness_scores = -np.sum(np.abs(new_population), axis=1)
# 
#                 # Calculate reward as improvement in average fitness (for minimization problem)
#                 reward = np.mean(self.fitness_scores) - np.mean(new_fitness_scores)
#                 reward = torch.tensor([reward], dtype=torch.float32)
# 
#                 # Update state and fitness
#                 self.population = new_population
#                 self.fitness_scores = new_fitness_scores
# 
#                 # Update state with new statistics
#                 pop_mean = self.population.mean(axis=0)
#                 pop_var = self.population.var(axis=0)
#                 best_ind = self.population[np.argmin(self.fitness_scores)]
#                 state_components = [pop_mean, pop_var, best_ind]
#                 state_vector = np.concatenate(state_components)
#                 next_state = torch.tensor(state_vector, dtype=torch.float32).unsqueeze(0)
# 
#                 # Store trajectory information
#                 log_probs.append(dist.log_prob(action))
#                 values.append(value)
#                 rewards.append(reward)
#                 entropies.append(entropy)
# 
#                 self.state = next_state
#                 self.episode_reward += reward.item()
# 
#                 self.global_step_counter.increment()
# 
#                 # Check if episode is done
#                 if self.global_step_counter.get() % MAX_EPISODE_LENGTH == 0:
#                     self.done = True
#                     break
# 
#             # If episode is done, use 0 as final state value
#             if self.done:
#                 R = torch.zeros(1, 1)
#                 self.episode_rewards.put(self.episode_reward)
#             else:
#                 # Bootstrap with value function
#                 _, R, _ = self.local_model(self.state)
# 
#             # Calculate returns and advantages
#             values.append(R)
#             R = R.detach()
# 
#             returns = []
#             for i in reversed(range(len(rewards))):
#                 R = rewards[i] + GAMMA * R
#                 returns.insert(0, R)
# 
#             returns = torch.cat(returns)
#             log_probs = torch.cat(log_probs)
#             values = torch.cat(values[:-1])
#             entropies = torch.stack(entropies)
# 
#             advantages = returns - values
# 
#             # Calculate losses
#             actor_loss = -(log_probs * advantages.detach()).mean()
#             critic_loss = F.mse_loss(values, returns)
#             entropy_loss = -entropies.mean()
# 
#             total_loss = actor_loss + 0.5 * critic_loss + ENTROPY_BETA * entropy_loss
# 
#             # Update global model
#             self.optimizer.zero_grad()
#             total_loss.backward()
# 
#             # Ensure gradient is in reasonable range
#             for param in self.local_model.parameters():
#                 if param.grad is not None:
#                     param.grad.data.clamp_(-1, 1)
# 
#             # Push local gradients to global model
#             for local_param, global_param in zip(self.local_model.parameters(), self.global_model.parameters()):
#                 if global_param.grad is None:
#                     global_param.grad = local_param.grad
#                 else:
#                     global_param.grad += local_param.grad
# 
#             self.optimizer.step()
# 
# # Communication queues
# fitness_queue = mp.Queue()
# result_queue = mp.Queue()
# 
# # Initialize and train the A3C model
# # Using expanded state size: mean + variance + best individual
# state_size = IND_SIZE * 3  # mean, variance, best individual
# global_model = ActorCriticWithAttention(state_size, POP_SIZE)
# global_model.share_memory()  # Share model parameters across processes
# optimizer = optim.Adam(global_model.parameters(), lr=LR)
# 
# # Initialize workers
# workers = [Worker(global_model, optimizer, global_step_counter, episode_rewards, i,
#                  fitness_queue, result_queue) 
#            for i in range(NUM_PROCESSES)]
# 
# # State for main process
# mean_rewards = deque(maxlen=100)
# training_complete = False
# best_reward = -float('inf')
# workers_started = False
# 
# # Initialize elite memory
# elite_memory = None
# elite_fitness = float('inf')
# 
# # Function to modify population using A3C and Attention
# def get_attention_modified_population(mat_population, fitness_values, bounds):
#     global workers, training_complete, mean_rewards, best_reward, global_model
#     global workers_started, elite_memory, elite_fitness, POP_SIZE, IND_SIZE, state_size
# 
#     # Convert input to proper numpy array
#     population = np.array(mat_population)
# 
#     # Handle the case where input is 1D (single individual)
#     if population.ndim == 1:
#         population = population.reshape(1, -1)
# 
#     # Get the actual dimensions
#     pop_size, ind_size = population.shape
# 
#     # Update global variables if dimensions change
#     if ind_size != IND_SIZE or pop_size != POP_SIZE:
#         POP_SIZE = pop_size
#         IND_SIZE = ind_size
#         state_size = ind_size * 3  # Update state size for new dimensions
# 
#     # Ensure dimensions match expected values
#     if ind_size != IND_SIZE or pop_size != POP_SIZE:
#         # Reinitialize model with new dimensions
#         global_model = ActorCriticWithAttention(state_size, pop_size)
#         global_model.share_memory()
#         optimizer = optim.Adam(global_model.parameters(), lr=LR)
# 
#         # Stop existing workers if any
#         if workers_started:
#             for worker in workers:
#                 worker.terminate()
#                 worker.join(timeout=0.1)
# 
#         # Create new workers
#         workers = [Worker(global_model, optimizer, global_step_counter, episode_rewards, i,
#                          fitness_queue, result_queue) 
#                   for i in range(NUM_PROCESSES)]
#         for worker in workers:
#             worker.daemon = True
# 
#         # Reset training state
#         training_complete = False
#         mean_rewards = deque(maxlen=100)
#         best_reward = -float('inf')
#         workers_started = False
# 
#         # Reset elite memory
#         elite_memory = None
#         elite_fitness = float('inf')
# 
#     # Start workers if not already started
#     if not workers_started:
#         for worker in workers:
#             worker.start()
#         workers_started = True
# 
#     # If fitness values are provided, put them in the queue
#     if fitness_values is not None:
#         try:
#             # Convert to numpy array if needed
#             fitness_values = np.array(fitness_values)
#             fitness_queue.put(fitness_values, block=False)
#         except:
#             pass  # Queue might be full, skip
# 
#     # Update elite memory if better solution found
#     if fitness_values is not None:
#         # Convert to numpy array if needed
#         fitness_values = np.array(fitness_values)
#         min_idx = np.argmin(fitness_values)
#         if fitness_values[min_idx] < elite_fitness:
#             elite_fitness = fitness_values[min_idx]
#             elite_memory = population[min_idx].copy()
# 
#     # Collect episode rewards if available
#     while not episode_rewards.empty():
#         try:
#             reward = episode_rewards.get(block=False)
#             mean_rewards.append(reward)
# 
#             # Check if we've reached good performance
#             if len(mean_rewards) >= 100:
#                 avg_reward = sum(mean_rewards) / len(mean_rewards)
#                 if avg_reward > best_reward:
#                     best_reward = avg_reward
#                     torch.save(global_model.state_dict(), 'best_a3c_model.pth')
# 
#                 if avg_reward > 10:  # Threshold for good performance
#                     training_complete = True
#         except:
#             break
# 
#     # Check if there's a result from workers
#     try:
#         new_population = result_queue.get(block=False)
#     except:
#         # If no result available, use the model to generate new population
#         with torch.no_grad():
#             # Create enhanced state representation
#             pop_mean = population.mean(axis=0)
#             pop_var = population.var(axis=0)
# 
#             if fitness_values is not None:
#                 fitness_values = np.array(fitness_values)
#                 best_idx = np.argmin(fitness_values)
#                 best_ind = population[best_idx]
#             else:
#                 # If no fitness provided, use L1 norm as proxy
#                 proxy_fitness = np.sum(np.abs(population), axis=1)
#                 best_idx = np.argmin(proxy_fitness)
#                 best_ind = population[best_idx]
# 
#             # Concatenate for rich state representation
#             state_components = [pop_mean, pop_var, best_ind]
#             state_vector = np.concatenate(state_components)
# 
#             state = torch.tensor(state_vector, dtype=torch.float32).unsqueeze(0)
#             action_probs, _, attention_weights = global_model(state)
#             attention_weights = attention_weights.numpy().flatten()
# 
#         # Ensure valid probability distribution
#         if np.isnan(attention_weights).any() or np.isinf(attention_weights).any():
#             attention_weights = np.ones(pop_size) / pop_size
#         else:
#             attention_weights = attention_weights / attention_weights.sum()
# 
#         # Generate new population through crossover and mutation
#         new_population = np.zeros_like(population)
#         for i in range(pop_size):
#             # Select parents using attention weights
#             parent_indices = np.random.choice(pop_size, 2, p=attention_weights)
#             parent1, parent2 = population[parent_indices[0]], population[parent_indices[1]]
# 
#             # Crossover
#             crossover_point = np.random.randint(1, ind_size)
#             child = np.concatenate([parent1[:crossover_point], parent2[crossover_point:]])
# 
#             # Adaptive mutation
#             mutation_prob = 0.1 if not training_complete else 0.05  # Reduce mutation when trained
#             mutation_mask = np.random.random(ind_size) < mutation_prob
#             child[mutation_mask] += np.random.normal(0, 0.1, size=mutation_mask.sum())
# 
#             new_population[i] = child
# 
#     # Include elite solution in new population
#     if elite_memory is not None:
#         # Replace worst individual with elite if we have fitness values
#         if fitness_values is not None:
#             fitness_values = np.array(fitness_values)
#             worst_idx = np.argmax(fitness_values)
#             new_population[worst_idx] = elite_memory.copy()
# 
#     # Apply boundary constraints if provided
#     if bounds is not None:
#         # Handle different types of bounds input
#         if isinstance(bounds, (list, tuple)) and len(bounds) == 2:
#             lb, ub = bounds
# 
#             # Convert to numpy arrays if they aren't already
#             if not isinstance(lb, np.ndarray):
#                 lb = np.array(lb)
#             if not isinstance(ub, np.ndarray):
#                 ub = np.array(ub)
# 
#             # Make sure they're the right shape
#             if lb.size == 1:
#                 lb = np.ones(ind_size) * lb
#             if ub.size == 1:
#                 ub = np.ones(ind_size) * ub
# 
#             # Reflective boundary handling
#             for i in range(pop_size):
#                 for j in range(ind_size):
#                     while new_population[i, j] < lb[j] or new_population[i, j] > ub[j]:
#                         if new_population[i, j] < lb[j]:
#                             new_population[i, j] = 2 * lb[j] - new_population[i, j]
#                         if new_population[i, j] > ub[j]:
#                             new_population[i, j] = 2 * ub[j] - new_population[i, j]
# 
#     return new_population.tolist()



from cec2017.functions import f1  # Import the fitness function f1 from cec2017
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
import torch.multiprocessing as mp
import numpy as np
import os
import time
from collections import deque

# Set up multiprocessing support
try:
    mp.set_start_method('spawn', force=True)
except RuntimeError:
    pass  # Already set

# Define Population & Attention Parameters
POP_SIZE = 30
IND_SIZE = 30
GAMMA = 0.99
LR = 0.0003
NUM_PROCESSES = os.cpu_count() - 1 or 1
MAX_EPISODE_LENGTH = 100
ENTROPY_BETA = 0.01

# Shared counter for tracking global steps
class Counter:
    def __init__(self):
        self.val = mp.Value('i', 0)
    
    def increment(self):
        with self.val.get_lock():
            self.val.value += 1
    
    def get(self):
        with self.val.get_lock():
            return self.val.value

global_step_counter = Counter()
episode_rewards = mp.Queue()

# Actor-Critic Network with Attention and Layer Normalization
class ActorCriticWithAttention(nn.Module):
    def __init__(self, state_dim, action_dim, hidden_dim=128):
        super(ActorCriticWithAttention, self).__init__()
        
        # Shared feature extraction layers
        self.shared_fc1 = nn.Linear(state_dim, hidden_dim)
        self.shared_fc2 = nn.Linear(hidden_dim, hidden_dim)
        
        # Layer normalization for stability
        self.layer_norm1 = nn.LayerNorm(hidden_dim)
        self.layer_norm2 = nn.LayerNorm(hidden_dim)
        
        # Self-attention mechanism
        self.query = nn.Linear(hidden_dim, hidden_dim)
        self.key = nn.Linear(hidden_dim, hidden_dim)
        self.value = nn.Linear(hidden_dim, hidden_dim)
        
        # Layer normalization for attention
        self.attention_norm = nn.LayerNorm(hidden_dim)
        
        # Actor (policy) head
        self.actor = nn.Linear(hidden_dim, action_dim)
        
        # Critic (value) head
        self.critic = nn.Linear(hidden_dim, 1)
    
    def forward(self, x):
        # Shared feature extraction with layer normalization
        x = F.relu(self.layer_norm1(self.shared_fc1(x)))
        x = F.relu(self.layer_norm2(self.shared_fc2(x)))
        
        # Self-attention
        batch_size = x.size(0)
        query = self.query(x)
        key = self.key(x)
        value = self.value(x)
        
        # Compute attention scores with temperature scaling for stability
        temperature = np.sqrt(key.size(1))
        if batch_size > 1:
            attention_scores = torch.bmm(query.unsqueeze(1), key.unsqueeze(2)).squeeze(1) / temperature
            attention_probs = F.softmax(attention_scores, dim=1)
            context = attention_probs.unsqueeze(1) * value
            context = context.sum(dim=1)
        else:
            attention_scores = torch.mm(query, key.t()) / temperature
            attention_probs = F.softmax(attention_scores, dim=1)
            context = torch.mm(attention_probs, value)
        
        # Combine with residual connection and normalization
        x = self.attention_norm(x + context)
        
        # Actor output (action probabilities)
        logits = self.actor(x)
        action_probs = F.softmax(logits, dim=-1)
        
        # Critic output (state value)
        state_value = self.critic(x)
        
        return action_probs, state_value, attention_probs

# Worker process for A3C
class Worker(mp.Process):
    def __init__(self, global_model, optimizer, global_step_counter, episode_rewards, 
                 process_idx, fitness_queue, result_queue):
        super(Worker, self).__init__()
        self.process_idx = process_idx
        self.global_model = global_model
        self.optimizer = optimizer
        self.global_step_counter = global_step_counter
        self.episode_rewards = episode_rewards
        self.fitness_queue = fitness_queue  # Queue to receive external fitness values
        self.result_queue = result_queue    # Queue to send back modified populations
        
        # Local model for this worker
        self.local_model = ActorCriticWithAttention(IND_SIZE, POP_SIZE)
        self.local_model.load_state_dict(self.global_model.state_dict())
        
        # Episode information
        self.done = True
        self.episode_reward = 0
        self.population = None
        self.state = None
        self.fitness_scores = None
        
    def run(self):
        while self.global_step_counter.get() < 10000:  # Limit total steps
            # Initialize or reset environment
            if self.done:
                # Generate new random population as initial state
                self.population = np.random.uniform(-100, 100, (POP_SIZE, IND_SIZE))
                
                # Use f1 as the fitness function (CEC 2017 benchmark)
                self.fitness_scores = f1(self.population)
                
                # Enhanced state representation: mean, variance, and best individual
                pop_mean = self.population.mean(axis=0)
                pop_var = self.population.var(axis=0)
                best_ind = self.population[np.argmin(self.fitness_scores)]
                
                # Concatenate statistics for richer state representation
                state_components = [pop_mean, pop_var, best_ind]
                state_vector = np.concatenate(state_components)
                
                self.state = torch.tensor(state_vector, dtype=torch.float32).unsqueeze(0)
                self.done = False
                self.episode_reward = 0
            
            self.local_model.load_state_dict(self.global_model.state_dict())
            
            # Storage for trajectory
            log_probs = []
            values = []
            rewards = []
            entropies = []
            
            # Collect trajectory
            for _ in range(MAX_EPISODE_LENGTH):
                # Get action probabilities and state value
                action_probs, value, attention_weights = self.local_model(self.state)
                
                # Sample action
                dist = torch.distributions.Categorical(action_probs)
                action = dist.sample()
                
                # Calculate entropy for exploration encouragement
                entropy = -(action_probs * action_probs.log()).sum()
                
                # Convert attention weights to numpy and ensure they're valid
                attention_weights = attention_weights.detach().numpy().flatten()
                if np.isnan(attention_weights).any() or np.isinf(attention_weights).any():
                    attention_weights = np.ones(POP_SIZE) / POP_SIZE
                
                # Generate new population through crossover and mutation
                new_population = np.zeros_like(self.population)
                for i in range(POP_SIZE):
                    # Select parents using attention weights
                    parent_indices = np.random.choice(POP_SIZE, 2, p=attention_weights/attention_weights.sum())
                    parent1, parent2 = self.population[parent_indices[0]], self.population[parent_indices[1]]
                    
                    # Crossover
                    crossover_point = np.random.randint(1, IND_SIZE)
                    child = np.concatenate([parent1[:crossover_point], parent2[crossover_point:]])
                    
                    # Adaptive mutation based on training progress
                    mutation_rate = max(0.05, 0.2 * (1 - self.global_step_counter.get() / 10000))
                    mutation_mask = np.random.random(IND_SIZE) < mutation_rate
                    mutation_strength = 0.1 * (1 - self.global_step_counter.get() / 10000)
                    child[mutation_mask] += np.random.normal(0, mutation_strength, size=mutation_mask.sum())
                    
                    new_population[i] = child
                
                # Use f1 as the fitness function (CEC 2017 benchmark)
                new_fitness_scores = f1(new_population)
                
                # Calculate reward as improvement in average fitness (for minimization problem)
                reward = np.mean(self.fitness_scores) - np.mean(new_fitness_scores)
                reward = torch.tensor([reward], dtype=torch.float32)
                
                # Update state and fitness
                self.population = new_population
                self.fitness_scores = new_fitness_scores
                
                # Update state with new statistics
                pop_mean = self.population.mean(axis=0)
                pop_var = self.population.var(axis=0)
                best_ind = self.population[np.argmin(self.fitness_scores)]
                state_components = [pop_mean, pop_var, best_ind]
                state_vector = np.concatenate(state_components)
                next_state = torch.tensor(state_vector, dtype=torch.float32).unsqueeze(0)
                
                # Store trajectory information
                log_probs.append(dist.log_prob(action))
                values.append(value)
                rewards.append(reward)
                entropies.append(entropy)
                
                self.state = next_state
                self.episode_reward += reward.item()
                
                self.global_step_counter.increment()
                
                # Check if episode is done
                if self.global_step_counter.get() % MAX_EPISODE_LENGTH == 0:
                    self.done = True
                    break
            
            # If episode is done, use 0 as final state value
            if self.done:
                R = torch.zeros(1, 1)
                self.episode_rewards.put(self.episode_reward)
