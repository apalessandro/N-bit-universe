# Visio

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

You can supply an arbitrary partition of the microscopic state space (all 2^N bitstrings) and treat its blocks as macrostates.

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
- `--layout`: Layout algorithm for the trajectory coarse graph when `--plot` is given. Choices: `spring` (default), `circular`, `shell`, `kamada`, `spectral`, `planar` (falls back to spring if not planar).
- `--coarse-node-size`: Node size for coarse-grained graph nodes when `--plot` is used (default 1200).

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
- `--layout`: Graph layout algorithm (default: `spring`). Other choices: `circular`, `shell`, `kamada`, `spectral`, `planar`.
- `--coarse-node-size`: Uniform node size for macrostates (default 1200) when using `coarse-graph` or `demo --plot`.

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
- `--layout`: Graph layout algorithm (default: `spring`). Other choices: `circular`, `shell`, `kamada`, `spectral`, `planar`.

Features:

- **Node labels**: Show macrostate value and microstate count
- **Edges**: Show transitions between macrostates
- **Edge labels**: Display the number of microstates making each transition (acts as the weight)
- Automatically tests for Markovian closure

### Graph Layouts

You can choose different layouts for visualization to emphasize structure:

- `spring`: Force-directed layout (Fruchterman-Reingold). Good general-purpose default.
- `circular`: Places nodes on a circle; useful for seeing cycle structure.
- `shell`: Concentric shells; NetworkX chooses grouping heuristics.
- `kamada`: Kamada–Kawai spring layout; sometimes produces more uniform edge lengths.
- `spectral`: Uses eigenvectors of graph Laplacian; can highlight connectivity clusters.
- `planar`: Attempts planar embedding (only works if graph is planar; falls back to `spring` otherwise).

If an unsupported layout string is provided or a planar embedding fails, the code automatically falls back to the `spring` layout.

## License

This project is provided as-is for educational and research purposes.
