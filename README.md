# N-bit Universe

A finite descriptive universe simulator implementing elementary cellular automata on N bits with comprehensive tools for analyzing microscopic and macroscopic dynamics.

## Overview

This project implements **any elementary cellular automaton** (rules 0-255) as finite, deterministic universes where:

- States are N-bit binary strings with periodic boundary conditions
- Evolution is synchronous and deterministic
- Multiple rule choices enable exploration of different dynamical behaviors
- The system can be analyzed at both microscopic (individual bits) and macroscopic (coarse-grained) levels

## Features

### Core Dynamics

- **Universal Rule Support**: Any elementary CA rule from 0 to 255 using Wolfram's encoding
- **Synchronous Evolution**: Deterministic cellular automaton updates with periodic boundaries
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
python n-bit_universe.py demo -N 4 -t 16 -c parity -r 90
```

Options:

- `-N`: Number of bits (default: 4)
- `-t, --steps`: Number of time steps (default: 16)
- `-r, --rule`: Cellular automaton rule number from **0 to 255** (default: 90)
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
N=4, steps=16, rule=90, coarse=parity
First 10 macrostates along trajectory: [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
Histogram of macrostates on trajectory: {1: 9, 0: 8}
Macro entropy over trajectory bag (nats): 0.692419
Markovianly closed? True
```

#### Cycle Analysis

Decompose the permutation into cycles:

```bash
python n-bit_universe.py cycles -N 4 -r 90
```

Options:

- `-N`: Number of bits (default: 4)
- `-r, --rule`: Cellular automaton rule number from **0 to 255** (default: 90)

#### Microscopic Phase Space Visualization

Visualize the complete state space graph showing all microscopic transitions:

```bash
python n-bit_universe.py graph -N 4 -r 90
```

Options:

- `-N`: Number of bits (default: 4)
- `-r, --rule`: Cellular automaton rule number from **0 to 255** (default: 90)

Shows a directed graph where:

- Nodes represent individual microstates (bit strings)
- Edges show the deterministic evolution under the selected rule
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
- **Edges**: Show transitions between macrostates
- **Edge labels**: Display the number of microstates making each transition (acts as the weight)
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
from n_bit_universe import (
    make_rule_step,  # General rule generator
    rule30_step,
    rule90_step,
    rule110_step,
    rule184_step,
    get_rule,
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

# Method 1: Use the general rule generator for any rule 0-255
rule_fn = make_rule_step(110)  # Creates rule-110

# Method 2: Use get_rule() for string-based lookup
rule_fn = get_rule("110")  # Also works with string

# Method 3: Use preset functions directly
rule_fn = rule110_step

# Evolve the system
trajectory = evolve(s0, steps=10, step_fn=rule_fn)

# Apply coarse-graining
macro_states = [cg_parity(s) for s in trajectory]

# Check if coarse-graining is Markovian
is_closed = is_markovianly_closed(N, cg_parity, rule_fn)

# Get induced macro transition matrix
macro_labels, M = induced_macro_transition_uniform(N, cg_parity, rule_fn)

# Analyze cycle structure
cycle_list = cycles(N, rule_fn)

# Explore any rule!
for rule_num in [22, 54, 60, 102, 150, 182]:
    rule = make_rule_step(rule_num)
    cycs = cycles(N, rule)
    print(f"Rule-{rule_num}: {len(cycs)} cycles")
```

## Theory

### Wolfram's Elementary Cellular Automaton Encoding

Every elementary CA rule is identified by a number from **0 to 255**. The rule number encodes the output for all 8 possible 3-bit neighborhoods:

```text
Neighborhood:  111  110  101  100  011  010  001  000
Bit position:   7    6    5    4    3    2    1    0
Output:      (rule_number >> bit_position) & 1
```

For example, **Rule-110** = 110₁₀ = 01101110₂:

```text
Neighborhood:  111  110  101  100  011  010  001  000
Output:         0    1    1    0    1    1    1    0
```

So when the neighborhood is `011`, the center cell becomes `1` (bit 3 of 110).

### Notable Rules

#### Rule-30 (Chaotic)

```text
s'[i] = s[i] ⊕ (s[i-1] OR s[i+1])
```

Exhibits chaotic behavior and is used in random number generation. Shows complex, unpredictable patterns.

#### Rule-90 (Additive)

```text
s'[i] = s[i-1] ⊕ s[i+1]
```

Produces patterns similar to Pascal's triangle modulo 2. Exhibits regular, symmetric structures. This rule is **additive** and has well-understood mathematical properties.

#### Rule-110 (Universal)

```text
s'[i] = s[i+1] ⊕ (s[i-1] OR s[i])
```

Proven to be **Turing complete** (capable of universal computation). Shows complex, structured behavior at the "edge of chaos."

#### Rule-184 (Traffic Flow)

```text
s'[i] = (s[i-1] AND s[i]) OR (s[i] AND NOT s[i+1]) OR (s[i-1] AND NOT s[i+1])
```

Models **traffic flow** where 1-bits represent cars. Cars move right if the space ahead is empty. Used in studying traffic dynamics and particle systems.

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

### Example 1: Exploring Arbitrary Rules

```python
from n_bit_universe import make_rule_step, evolve, int_to_state, cycles

N = 6
s0 = int_to_state(1, N)

# Systematically explore rules
for rule_num in range(0, 256, 10):  # Sample every 10th rule
    rule_fn = make_rule_step(rule_num)
    cycs = cycles(N, rule_fn)
    lengths = sorted(len(c) for c in cycs)
    print(f"Rule-{rule_num:3d}: {len(cycs):2d} cycles, lengths {lengths[:5]}...")
```

### Example 2: Comparing Rules

```python
from n_bit_universe import get_rule, evolve, int_to_state, cycles

N = 6
s0 = int_to_state(1, N)

# Compare cycle structures across rules
for rule_name in ["30", "90", "110", "184"]:
    rule_fn = get_rule(rule_name)
    cycs = cycles(N, rule_fn)
    lengths = sorted(len(c) for c in cycs)
    print(f"Rule-{rule_name}: {len(cycs)} cycles with lengths {lengths}")
```

### Example 3: Checking Different Coarse-Grainings

```python
from n_bit_universe import (
    cg_parity, cg_weight, cg_rotation_class,
    is_markovianly_closed, get_rule
)

N = 5
rule_fn = get_rule("110")

# Check which coarse-grainings are Markovian for rule-110
print(f"Parity closed: {is_markovianly_closed(N, cg_parity, rule_fn)}")
print(f"Weight closed: {is_markovianly_closed(N, cg_weight, rule_fn)}")
print(f"Rotation closed: {is_markovianly_closed(N, cg_rotation_class, rule_fn)}")
```

### Example 4: Exploring Trajectories

```python
from n_bit_universe import evolve, int_to_state, state_to_int, get_rule

# Compare trajectory periods for different rules
N = 6
rule_fn = get_rule("184")  # Traffic flow rule

for initial_value in [1, 7, 15, 63]:
    s0 = int_to_state(initial_value, N)
    traj = evolve(s0, steps=50, step_fn=rule_fn)
    
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
