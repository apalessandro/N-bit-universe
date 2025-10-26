# Usage Examples for N-bit Universe

## General Rule System

This simulator now supports **all 256 elementary cellular automaton rules** (0-255) using Wolfram's encoding.

### Quick Start: Try Any Rule

```bash
# Rule-30: Chaotic
python n-bit_universe.py demo -N 5 -t 20 -r 30 -c parity

# Rule-90: Additive (default)
python n-bit_universe.py demo -N 5 -t 20 -r 90 -c parity

# Rule-110: Turing complete
python n-bit_universe.py demo -N 5 -t 20 -r 110 -c parity

# Rule-184: Traffic flow
python n-bit_universe.py demo -N 5 -t 20 -r 184 -c parity

# Explore lesser-known rules!
python n-bit_universe.py demo -N 5 -t 20 -r 22 -c parity
python n-bit_universe.py demo -N 5 -t 20 -r 54 -c parity
python n-bit_universe.py demo -N 5 -t 20 -r 150 -c parity
```

### Systematic Rule Exploration

```bash
# Compare cycle structures
python n-bit_universe.py cycles -N 6 -r 0
python n-bit_universe.py cycles -N 6 -r 50
python n-bit_universe.py cycles -N 6 -r 100
python n-bit_universe.py cycles -N 6 -r 150
python n-bit_universe.py cycles -N 6 -r 200
python n-bit_universe.py cycles -N 6 -r 255

# Visualize phase space for different rules
python n-bit_universe.py graph -N 4 -r 30
python n-bit_universe.py graph -N 4 -r 110
python n-bit_universe.py graph -N 4 -r 126
```

### Python API: Generate Any Rule

```python
from n_bit_universe import make_rule_step, evolve, int_to_state, cycles

N = 6

# Scan through all rules
for rule_num in range(256):
    rule_fn = make_rule_step(rule_num)
    cycs = cycles(N, rule_fn)
    
    # Find rules with interesting properties
    if len(cycs) == 1:  # Single cycle = permutation
        print(f"Rule-{rule_num}: Single cycle of length {len(cycs[0])}")
    elif len(cycs) > 20:  # Many small cycles
        print(f"Rule-{rule_num}: {len(cycs)} cycles")
```

### Understanding Wolfram Encoding

Rule numbers encode the lookup table for all 8 possible neighborhoods:

```python
def explain_rule(rule_num):
    """Show what a rule does for each neighborhood."""
    neighborhoods = [
        (1,1,1), (1,1,0), (1,0,1), (1,0,0),
        (0,1,1), (0,1,0), (0,0,1), (0,0,0)
    ]
    
    print(f"Rule-{rule_num}:")
    for i, (l, c, r) in enumerate(neighborhoods):
        output = (rule_num >> i) & 1
        print(f"  {l}{c}{r} → {output}")

# Example: Rule-110
explain_rule(110)
# Output:
# Rule-110:
#   111 → 0
#   110 → 1
#   101 → 1
#   100 → 0
#   011 → 1
#   010 → 1
#   001 → 1
#   000 → 0
```

### Finding Rules by Property

```python
from n_bit_universe import make_rule_step, is_markovianly_closed, cg_parity

N = 5

# Find rules where parity is Markovian
markov_rules = []
for rule_num in range(256):
    rule_fn = make_rule_step(rule_num)
    if is_markovianly_closed(N, cg_parity, rule_fn):
        markov_rules.append(rule_num)

print(f"Rules with Markovian parity closure: {markov_rules}")
```

### Comparing Rule Classes

Wolfram classified rules into 4 classes:

- **Class 1**: Converge to uniform state (e.g., 0, 32, 64)
- **Class 2**: Simple stable/periodic patterns (e.g., 4, 108, 170)
- **Class 3**: Chaotic/random (e.g., 30, 45, 106)
- **Class 4**: Complex, long-lived structures (e.g., 110, 124)

```bash
# Compare classes
python n-bit_universe.py demo -N 8 -t 30 -r 0   -c weight  # Class 1
python n-bit_universe.py demo -N 8 -t 30 -r 4   -c weight  # Class 2
python n-bit_universe.py demo -N 8 -t 30 -r 30  -c weight  # Class 3
python n-bit_universe.py demo -N 8 -t 30 -r 110 -c weight  # Class 4
```

## Advanced Usage

### Custom Analysis Script

```python
from n_bit_universe import make_rule_step, evolve, int_to_state, cycles
import numpy as np

def analyze_rule(rule_num, N=6):
    """Comprehensive analysis of a rule."""
    rule_fn = make_rule_step(rule_num)
    
    # Cycle structure
    cycs = cycles(N, rule_fn)
    max_cycle = max(len(c) for c in cycs)
    
    # Sample trajectory
    s0 = int_to_state(1, N)
    traj = evolve(s0, steps=50, step_fn=rule_fn)
    
    # Check for fixed points
    fixed_points = sum(1 for c in cycs if len(c) == 1)
    
    return {
        'rule': rule_num,
        'num_cycles': len(cycs),
        'max_cycle_length': max_cycle,
        'fixed_points': fixed_points,
        'trajectory_sample': traj[:5]
    }

# Analyze interesting rules
for rule in [30, 90, 110, 184]:
    result = analyze_rule(rule)
    print(f"\nRule-{rule}:")
    print(f"  Cycles: {result['num_cycles']}")
    print(f"  Max cycle length: {result['max_cycle_length']}")
    print(f"  Fixed points: {result['fixed_points']}")
```

### Visualize Rule Space

```bash
# Create a comparative visualization
for r in 22 54 60 90 110 126 150 182; do
    python n-bit_universe.py coarse-graph -N 6 -c parity -r $r
done
```

### NEW: Custom Coarse-Graining Examples

Inline JSON specifying a partition (N=3 has 8 states 000..111):

```bash
python n-bit_universe.py coarse-graph -N 3 -c custom --groups '{"Group1":["000","001"],"Group2":["010","011"],"Rest":["100","101","110","111"]}'
```

Using integer encodings mixed with bitstrings:

```bash
python n-bit_universe.py demo -N 3 -t 10 -c custom --groups '{"A":[0,1],"B":["010","011"],"C":[4,5,6,7]}'
```

List form (labels auto-assigned as 0,1,2,...):

```bash
python n-bit_universe.py demo -N 3 -t 10 -c custom --groups '[["000","001"],["010","011"],["100","101","110","111"]]'
```

From file:

```bash
python n-bit_universe.py coarse-graph -N 4 -c custom --groups custom_partition_example.json
```

If a state is missing or duplicated you'll get an error pointing out missing / duplicate assignments so you can fix the partition quickly.
