# N-bit Universe

A finite descriptive universe simulator implementing the rule-90 cellular automaton on N bits with comprehensive tools for analyzing microscopic and macroscopic dynamics.

## Overview

This project implements a **rule-90 cellular automaton** as a finite, deterministic universe where:

- States are N-bit binary strings with periodic boundary conditions
- Evolution is synchronous and deterministic
- The system can be analyzed at both microscopic (individual bits) and macroscopic (coarse-grained) levels

## Features

### Core Dynamics

- **Rule-90 Evolution**: Synchronous cellular automaton updates with periodic boundaries
- **Trajectory Generation**: Track system evolution over time
- **Permutation Analysis**: View the state space as a permutation group with cycle decomposition

### Coarse-Graining Methods

Three built-in macroscopic descriptions:

1. **Hamming Weight**: Count of 1-bits in the state
2. **Parity**: Even/odd number of 1-bits
3. **Rotation Class**: Equivalence class under cyclic rotations

### Analysis Tools

- **Entropy Calculations**: Shannon entropy for microscopic and macroscopic distributions (natural logs)
- **Markovian Closure Test**: Check if a coarse-graining preserves Markov structure (lumpability)
- **Induced Macro Transitions**: Build transition operators at the macroscopic level
- **Cycle Decomposition**: Analyze the permutation structure of the dynamics

## Installation

Requires Python 3.7+ and NumPy:

```bash
pip install numpy
```

## Usage

### Command Line Interface

#### Demo Mode

Run a simulation with macro-level analysis and visualization:

```bash
python n-bit_universe.py demo -N 4 -t 16 -c parity
```

Options:

- `-N`: Number of bits (default: 4)
- `-t, --steps`: Number of time steps (default: 16)
- `-c, --coarse`: Coarse-graining method: `parity`, `weight`, or `rotation` (default: parity)

The demo mode provides:

- **Trajectory visualization**: Shows the path through the coarse-grained state space
  - Green node marks the starting macrostate
  - Red edges highlight the trajectory path
  - Dark red nodes show visited macrostates
- **Entropy tracking**: Plots the macrostate entropy at each time step
- **Markovian closure test**: Checks if the coarse-graining preserves Markov structure
- **Transition matrix**: Shows induced macro-level dynamics

Example output:

```text
N=4, steps=16, coarse=parity
First 10 macrostates along trajectory: [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
Histogram of macrostates on trajectory: {1: 9, 0: 8}
Macro entropy over trajectory bag (nats): 0.692419
Markovianly closed? True
```

#### Cycle Analysis

Decompose the permutation into cycles:

```bash
python n-bit_universe.py cycles -N 4
```

Options:

- `-N`: Number of bits (default: 4)

#### Microscopic Phase Space Visualization

Visualize the complete state space graph showing all microscopic transitions:

```bash
python n-bit_universe.py graph -N 4
```

Options:

- `-N`: Number of bits (default: 4)

Shows a directed graph where:

- Nodes represent individual microstates (bit strings)
- Edges show the deterministic evolution under rule-90
- Cycles are visible in the graph structure

#### Coarse-Grained Phase Space Visualization

Visualize the macroscopic state space with different coarse-graining methods:

```bash
python n-bit_universe.py coarse-graph -N 8 -c rotation
```

Options:

- `-N`: Number of bits (default: 4)
- `-c, --coarse`: Coarse-graining method: `parity`, `weight`, or `rotation` (default: parity)

Features:

- **Node size**: Proportional to the number of microstates in each macrostate (capped at maximum size)
- **Node labels**: Show macrostate value and microstate count
- **Edge width**: Proportional to transition probability
- **Edge labels**: Display the number of microstates making each transition
- Automatically tests for Markovian closure

Example:

```bash
# Visualize with parity coarse-graining
python n-bit_universe.py coarse-graph -N 16 -c parity

# Visualize with Hamming weight
python n-bit_universe.py coarse-graph -N 8 -c weight

# Visualize with rotation classes
python n-bit_universe.py coarse-graph -N 8 -c rotation
```

### Python API

```python
import numpy as np
from n-bit_universe import (
    rule90_step,
    evolve,
    cg_weight,
    cg_parity,
    cg_rotation_class,
    entropy_from_counts,
    is_markovianly_closed,
    induced_macro_transition_uniform,
    cycles,
    int_to_state,
    state_to_int
)

# Create initial state (e.g., 0001)
N = 4
s0 = int_to_state(1, N)

# Evolve the system
trajectory = evolve(s0, steps=10, step_fn=rule90_step)

# Apply coarse-graining
macro_states = [cg_parity(s) for s in trajectory]

# Check if coarse-graining is Markovian
is_closed = is_markovianly_closed(N, cg_parity, rule90_step)

# Get induced macro transition matrix
macro_labels, M = induced_macro_transition_uniform(N, cg_parity, rule90_step)

# Analyze cycle structure
cycle_list = cycles(N, rule90_step)
```

## Theory

### Rule-90 Cellular Automaton

Each bit updates according to:

```text
s'[i] = s[i-1] ⊕ s[i+1]
```

where ⊕ is XOR and indices wrap around (periodic boundaries).

### Coarse-Graining and Emergence

This simulator explores how macroscopic descriptions relate to microscopic dynamics:

- **Markovian Closure**: A coarse-graining is "closed" if the macro dynamics form a proper Markov chain
- **Entropy**: Measures information content at different scales
- **Induced Transitions**: How microscopic dynamics translate to macro-level probabilities

### Permutation Structure

For finite N, the dynamics form a permutation of the 2^N possible states. The cycle decomposition reveals:

- Period of trajectories
- Symmetry structure
- Conservation properties

## Examples

### Example 1: Checking Different Coarse-Grainings

```python
N = 5

# Check if parity is Markovian
print(f"Parity closed: {is_markovianly_closed(N, cg_parity)}")

# Check if weight is Markovian
print(f"Weight closed: {is_markovianly_closed(N, cg_weight)}")

# Check if rotation class is Markovian
print(f"Rotation closed: {is_markovianly_closed(N, cg_rotation_class)}")
```

### Example 2: Exploring Trajectories

```python
# Start from different initial conditions
N = 6
for initial_value in [1, 7, 15, 63]:
    s0 = int_to_state(initial_value, N)
    traj = evolve(s0, steps=20)
    
    # Find period
    visited = {}
    for t, s in enumerate(traj):
        key = state_to_int(s)
        if key in visited:
            period = t - visited[key]
            print(f"Initial {initial_value}: period {period}")
            break
        visited[key] = t
```

## Mathematical Notation

- **N**: Number of bits
- **s**: Microscopic state (N-bit binary string)
- **H**: Shannon entropy (in natural logarithms/nats)
- **M[m',m]**: Induced macro transition matrix (column-stochastic)
- **T**: Time evolution operator (rule-90 step)

## License

This project is provided as-is for educational and research purposes.

## Contributing

Feel free to extend this simulator with:

- Additional cellular automaton rules
- New coarse-graining functions
- Different boundary conditions
- Visualization tools
- Performance optimizations for larger N

## References

Rule-90 is a well-studied elementary cellular automaton (Wolfram's classification). This implementation focuses on:

- Finite universes (small N)
- Information theory perspectives
- Coarse-graining and emergence
- Markovian closure (lumpability in stochastic processes)
