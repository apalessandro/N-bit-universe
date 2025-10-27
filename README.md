# N-bit Universe

A finite descriptive universe simulator implementing elementary cellular automata on N bits with comprehensive tools for analyzing microscopic and macroscopic dynamics.

## Overview

This project implements **any deterministic rule** on N-bit states as finite universes where:

- States are N-bit binary strings
- Evolution is synchronous and deterministic
- Rules are specified as mappings from each microstate to its image
- The system can be analyzed at both microscopic (individual bits) and macroscopic (coarse-grained) levels

## Features

### Core Dynamics

- **General Rule Support**: Any deterministic rule specified as space-separated integers 0…2^N−1 representing the mapping from each microstate to its image (e.g., "2 3 0 1" for N=2)
- **Elementary Cellular Automaton Support**: Built-in `--eca` flag to automatically generate the full 2^N mapping from any ECA rule number (0-255), eliminating the need to hand-type permutations for standard rules like Rule 90, Rule 110, etc.
- **Synchronous Evolution**: Deterministic cellular automaton updates with periodic boundaries
- **Trajectory Generation**: Track system evolution over time
- **Permutation Analysis**: View the state space as a permutation group with cycle decomposition

### Coarse-Graining

You can supply an arbitrary partition of the microscopic state space (all 2^N bitstrings) and treat its blocks as macrostates. This lets you explore descriptive phase spaces.

Add `--groups <SPEC>` where `<SPEC>` is either:

1. A path to a JSON file, or
2. An inline JSON string.

Accepted JSON formats:

**Mapping form** (explicit labels):

```json
{
    "Low": ["000", "001", 2],
    "Mid": ["011", "010"],
    "High": [4, 5, 6, 7]
}
```

**List form** (labels become "0", "1", ... automatically):

```json
[
    ["000", "001", "010"],
    ["011", "100"],
    ["101", "110", "111"]
]
```

Elements inside groups can be:

- Bit strings of length N (e.g. `"0101"` for N=4)
- Integers (e.g. `13`) or numeric strings (`"13"`) representing the integer encoding of the bitstring.

Validation rules:

- Every one of the `2^N` microstates must appear **exactly once** (full partition)
- Duplicates or omissions raise an error with helpful hints
- Labels are arbitrary strings (mapping form) or auto-generated indices (list form)

Example (inline JSON):

```bash
python n-bit_universe.py coarse-graph -N 3 --groups '{"A":["000","001"],"B":["010","011"],"C":["100","101","110","111"]}'
```

Example (file):

```bash
python n-bit_universe.py coarse-graph -N 4 --groups custom_partition_example.json
```

You can also use `demo` mode:

```bash
python n-bit_universe.py demo -N 3 -t 12 -r "0 1 2 3 4 5 6 7" --groups '{"even-weight":[0,3,5,6],"odd-weight":[1,2,4,7]}' --plot
```

Returned macro labels (e.g. `"A"`, `"High"`, or `"2"`) integrate seamlessly with entropy computations, transition matrix construction, Markovian closure tests, and visualization.

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
python n-bit_universe.py demo -N 4 -t 16 -r "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15" --groups custom_partition_example.json

# Or use the --eca flag for elementary cellular automaton rules:
python n-bit_universe.py demo -N 4 -t 16 --eca 90 --groups custom_partition_example.json
```

Options:

- `-N`: Number of bits (default: 4)
- `-t, --steps`: Number of time steps (default: 16)
- `-r, --rule`: Space-separated integers 0…2^N−1 representing the mapping from each microstate to its image (e.g. "2 3 0 1" for N=2)
- `--eca RULE_NUM`: Elementary cellular automaton rule number (0-255). Automatically generates the full 2^N mapping. Mutually exclusive with `-r`.
- `--groups`: JSON string or path specifying a full partition of the 2^N microstates (required)

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
N=4, steps=16, rule=2 3 0 1
First 10 macrostates along trajectory: [1, 0, 1, 0, 1, 0, 1, 0, 1, 0]
Histogram of macrostates on trajectory: {1: 9, 0: 8}
Macro entropy over trajectory bag (nats): 0.692419
Markovianly closed? True
```

#### Cycle Analysis

Decompose the permutation into cycles:

```bash
python n-bit_universe.py cycles -N 4 -r "0 1 2 3"

# Or use --eca for elementary cellular automaton rules:
python n-bit_universe.py cycles -N 4 --eca 90
```

Options:

- `-N`: Number of bits (default: 4)
- `-r, --rule`: Space-separated integers 0…2^N−1 representing the mapping from each microstate to its image
- `--eca RULE_NUM`: Elementary cellular automaton rule number (0-255). Mutually exclusive with `-r`.

#### Microscopic Phase Space Visualization

Visualize the complete state space graph showing all microscopic transitions:

```bash
python n-bit_universe.py graph -N 4 -r "0 1 2 3"

# Or use --eca for elementary cellular automaton rules:
python n-bit_universe.py graph -N 4 --eca 110
```

Options:

- `-N`: Number of bits (default: 4)
- `-r, --rule`: Space-separated integers 0…2^N−1 representing the mapping from each microstate to its image
- `--eca RULE_NUM`: Elementary cellular automaton rule number (0-255). Mutually exclusive with `-r`.

Shows a directed graph where:

- Nodes represent individual microstates (bit strings)
- Edges show the deterministic evolution under the selected rule
- Cycles are visible in the graph structure

#### Coarse-Grained Phase Space Visualization

Visualize the macroscopic state space with custom coarse-graining:

```bash
python n-bit_universe.py coarse-graph -N 8 -r "0 1 2 3 4 5 6 7" --groups custom_partition_example.json

# Or use --eca for elementary cellular automaton rules:
python n-bit_universe.py coarse-graph -N 8 --eca 30 --groups custom_partition_example.json
```

Options:

- `-N`: Number of bits (default: 4)
- `-r, --rule`: Space-separated integers 0…2^N−1 representing the mapping from each microstate to its image (e.g. "0 1 2 3 4 5 6 7" for N=3)
- `--eca RULE_NUM`: Elementary cellular automaton rule number (0-255). Mutually exclusive with `-r`.
- `--groups`: JSON string or path specifying a full partition (required)

Features:

- **Node size**: Proportional to the number of microstates in each macrostate (capped at maximum size)
- **Node labels**: Show macrostate value and microstate count
- **Edges**: Show transitions between macrostates
- **Edge labels**: Display the number of microstates making each transition (acts as the weight)
- Automatically tests for Markovian closure

Examples:

```bash
# Visualize with a custom partition (inline JSON)
python n-bit_universe.py coarse-graph -N 3 -r "0 1 2 3 4 5 6 7" --groups '{"left":["000","001"],"middle":["010","011"],"right":["100","101","110","111"]}'

# Visualize with a custom partition from file
python n-bit_universe.py coarse-graph -N 4 -r "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15" --groups custom_partition_example.json

# Using --eca with Rule 90
python n-bit_universe.py coarse-graph -N 5 --eca 90 --groups custom_partition_example.json
```

### Using the --eca Flag

The `--eca` flag provides a convenient way to use elementary cellular automaton rules without manually typing the full permutation string. Simply specify the ECA rule number (0-255) and the system automatically generates the complete 2^N→2^N mapping.

**Example - Rule 90:**

```bash
# Instead of manually typing the full permutation:
python n-bit_universe.py cycles -N 5 -r "0 15 30 17 60 51 34 49 120 119 102 113 68 83 98 85"

# Simply use:
python n-bit_universe.py cycles -N 5 --eca 90
```

The `--eca` flag works with all commands (`demo`, `cycles`, `graph`, `coarse-graph`) and is mutually exclusive with the `-r/--rule` option.

## Theory

### General Rule Encoding

Every deterministic rule is specified as space-separated integers 0…2^N−1 representing the mapping from each microstate to its image.

For example, for N=2, the string `"2 3 0 1"` means:

- 0 → 2
- 1 → 3
- 2 → 0
- 3 → 1

For N=4 (16 states), you would use: `"0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15"` for the identity permutation.

This allows you to explore any deterministic mapping, including reversible rules (permutations) and non-reversible ones.

### Notable Examples

#### Example Permutation (Reversible)

For N=3, the string `"0 1 2 3 4 5 6 7"` is the identity permutation (each state maps to itself).

For N=2, the string `"2 3 0 1"` is a permutation with two cycles: (0→2→0) and (1→3→1).

For N=4, the string `"15 14 13 12 11 10 9 8 7 6 5 4 3 2 1 0"` reverses the order of all states.

#### Example Non-Reversible Rule

For N=2, the string `"0 0 0 0"` maps every state to 0 (not reversible).

You can specify any mapping you like, as long as each value is a valid state index (0 to 2^N-1).

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
from n_bit_universe import parse_permutation_rule, evolve, int_to_state, cycles
import itertools

N = 2  # Small N for demonstration
s0 = int_to_state(1, N)

# Generate a few example permutations
M = 2 ** N
example_rules = ["0 1 2 3", "3 2 1 0", "1 0 3 2", "2 3 0 1"]  # Some permutations for N=2

for rule_str in example_rules:
    rule_fn = parse_permutation_rule(rule_str, N)
    cycs = cycles(N, rule_fn)
    lengths = sorted(len(c) for c in cycs)
    print(f"Rule '{rule_str}': {len(cycs)} cycles, lengths {lengths}")
```

### Example 2: Comparing Rules

```python
from n_bit_universe import parse_permutation_rule, evolve, int_to_state, cycles

N = 3
s0 = int_to_state(1, N)

# Compare cycle structures across different rules
rule_examples = {
    "identity": "0 1 2 3 4 5 6 7",
    "reverse": "7 6 5 4 3 2 1 0",
    "rotation": "1 2 3 4 5 6 7 0",
}

for name, rule_str in rule_examples.items():
    rule_fn = parse_permutation_rule(rule_str, N)
    cycs = cycles(N, rule_fn)
    lengths = sorted(len(c) for c in cycs)
    print(f"Rule '{name}': {len(cycs)} cycles with lengths {lengths}")
```

### Example 3: Checking Different Coarse-Grainings

```python
from n_bit_universe import (
    is_markovianly_closed, parse_permutation_rule, _parse_custom_partition
)

N = 5
rule_str = "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31"  # Identity for N=5
rule_fn = parse_permutation_rule(rule_str, N)

# Define different custom partitions
parity_partition = '{"even":[0,3,5,6,9,10,12,15,17,18,20,23,24,27,29,30],"odd":[1,2,4,7,8,11,13,14,16,19,21,22,25,26,28,31]}'
weight_partition = '{"w0":[0],"w1":[1,2,4,8,16],"w2":[3,5,6,9,10,12,17,18,20,24],"w3":[7,11,13,14,19,21,22,25,26,28],"w4":[15,23,27,29,30],"w5":[31]}'

cg_parity = _parse_custom_partition(parity_partition, N)
cg_weight = _parse_custom_partition(weight_partition, N)

# Check which coarse-grainings are Markovian
print(f"Parity closed: {is_markovianly_closed(N, cg_parity, rule_fn)}")
print(f"Weight closed: {is_markovianly_closed(N, cg_weight, rule_fn)}")
```

### Example 4: Exploring Trajectories

```python
from n_bit_universe import evolve, int_to_state, state_to_int, parse_permutation_rule

# Compare trajectory periods for different rules
N = 4
rule_str = "0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15"  # Identity permutation
rule_fn = parse_permutation_rule(rule_str, N)

for initial_value in [1, 7, 15, 10]:
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
- **T**: Time evolution operator

## License

This project is provided as-is for educational and research purposes.
