"""
Finite Descriptive Universe: elementary cellular automaton on N bits.
Implements:
  - Deterministic synchronous updates with periodic boundary conditions
  - Any elementary cellular automaton rule (0-255) via Wolfram encoding
  - Custom coarse-grainings via JSON specification
  - Microscopic and macroscopic entropy (natural logs)
  - Induced macro transition operator (uniform over microstates within a macrocell)
  - Test for Markovian closure (lumpability)
  - Cycle decomposition of the microscopic permutation
  - Simple CLI for quick experiments
"""

import argparse
import json
import math
import os
from collections import Counter, defaultdict
from typing import Callable, Dict, Iterable, List, Sequence, Tuple

import numpy as np
import matplotlib.pyplot as plt
import networkx as nx

Array = np.ndarray
State = Array  # shape (N,), dtype uint8 with entries in {0,1}
StepFn = Callable[[State], State]
CoarseFn = Callable[[State], object]

# ----------------------------
# Microscopic dynamics
# ----------------------------


def make_rule_step(rule_number: int) -> StepFn:
    """
    Create a step function for any elementary cellular automaton (0-255).
    Uses Wolfram's rule encoding where the rule number's binary representation
    gives the output for each of the 8 possible 3-bit neighborhoods.

    Neighborhood encoding (left, center, right):
      111 -> bit 7, 110 -> bit 6, 101 -> bit 5, 100 -> bit 4,
      011 -> bit 3, 010 -> bit 2, 001 -> bit 1, 000 -> bit 0

    Args:
        rule_number: Integer from 0 to 255

    Returns:
        Step function that applies the rule
    """
    if not 0 <= rule_number <= 255:
        raise ValueError(f"Rule number must be 0-255, got {rule_number}")

    # Precompute lookup table from rule number
    lookup = [(rule_number >> i) & 1 for i in range(8)]

    def step(s: State) -> State:
        left = np.roll(s, 1)
        center = s
        right = np.roll(s, -1)

        # Compute neighborhood index: 4*left + 2*center + right
        neighborhood = 4 * left + 2 * center + right

        # Apply lookup table
        result = np.array([lookup[idx] for idx in neighborhood], dtype=np.uint8)
        return result

    return step


def generate_eca_rule_string(rule_number: int, N: int) -> str:
    """
    Generate the full permutation string for an elementary cellular automaton rule.

    Applies the ECA rule to all 2^N possible states to build the complete mapping.
    This allows using ECA rules (0-255) directly without hand-typing the full permutation.

    Args:
        rule_number: ECA rule number (0-255)
        N: Number of bits

    Returns:
        Space-separated string of 2^N integers representing the complete mapping

    Example:
        For Rule 90, N=3: generates "0 7 6 1 4 3 2 5"
        (Each state maps to the result of applying Rule 90 once)
    """
    if not 0 <= rule_number <= 255:
        raise ValueError(f"ECA rule number must be 0-255, got {rule_number}")

    step_fn = make_rule_step(rule_number)
    mapping = []

    for i in range(2**N):
        s = int_to_state(i, N)
        s_next = step_fn(s)
        mapping.append(state_to_int(s_next))

    return " ".join(map(str, mapping))


def parse_permutation_rule(rule_str: str, N: int) -> StepFn:
    """
    Parse a user-supplied string representing the mapping from each microstate to its image.

    Accepts space-separated integers 0…2^N−1: "2 3 0 1" (for any N)

    Example: For N=2, rule_str="2 3 0 1" means:
      0 -> 2
      1 -> 3
      2 -> 0
      3 -> 1

    Returns a step function.
    """
    M = 2**N

    parts = rule_str.strip().split()
    if len(parts) != M:
        raise ValueError(
            f"Rule string must have {M} space-separated integers for N={N}, got {len(parts)}"
        )
    try:
        mapping = [int(p) for p in parts]
    except ValueError as e:
        raise ValueError(f"Invalid integer in rule string: {e}")

    if not all(0 <= v < M for v in mapping):
        raise ValueError(f"All values in rule string must be in 0..{M - 1}")

    def step(s: State) -> State:
        idx = state_to_int(s)
        next_idx = mapping[idx]
        return int_to_state(next_idx, N)

    return step


def evolve(s0: State, steps: int, step_fn: StepFn) -> List[State]:
    """Return the trajectory [s0, s1, ..., s_steps]."""
    traj = [s0.copy()]
    s = s0.copy()
    for _ in range(steps):
        s = step_fn(s)
        traj.append(s.copy())
    return traj


# ----------------------------
# Entropy utilities (natural logs)
# ----------------------------


def entropy_from_counts(counts: Dict[object, int]) -> float:
    """Shannon entropy H = -sum p log p with natural logs."""
    total = sum(counts.values())
    if total == 0:
        return 0.0
    H = 0.0
    for c in counts.values():
        if c > 0:
            p = c / total
            H -= p * math.log(p)
    return H


def entropy_from_probs(probs: Iterable[float]) -> float:
    H = 0.0
    for p in probs:
        if p > 0.0:
            H -= p * math.log(p)
    return H


# ----------------------------
# State space enumeration (small N)
# ----------------------------


def int_to_state(x: int, N: int) -> State:
    return np.array([(x >> i) & 1 for i in range(N)][::-1], dtype=np.uint8)


def state_to_int(s: State) -> int:
    out = 0
    for bit in s.astype(int):
        out = (out << 1) | bit
    return out


def enumerate_states(N: int) -> List[State]:
    return [int_to_state(x, N) for x in range(2**N)]


def permutation_graph(N: int, step_fn: StepFn) -> Dict[int, int]:
    """Return the permutation as a dict: x -> y in integer labels."""
    mapping: Dict[int, int] = {}
    for s in enumerate_states(N):
        x = state_to_int(s)
        y = state_to_int(step_fn(s))
        mapping[x] = y
    return mapping


def cycles(N: int, step_fn: StepFn) -> List[List[int]]:
    """Decompose the permutation into cycles, states labeled as integers."""
    mapping = permutation_graph(N, step_fn)
    seen = set()
    out: List[List[int]] = []
    for start in mapping:
        if start in seen:
            continue
        cyc = []
        x = start
        while x not in seen:
            seen.add(x)
            cyc.append(x)
            x = mapping[x]
        out.append(cyc)
    return out


# ----------------------------
# Macro distributions and transitions
# ----------------------------


def macro_histogram(traj: Sequence[State], coarse: CoarseFn) -> Dict[object, int]:
    return Counter([coarse(s) for s in traj])


def induced_macro_transition_uniform(
    N: int,
    coarse: CoarseFn,
    step_fn: StepFn,
) -> Tuple[List[object], np.ndarray]:
    """
    Build the induced macro transition operator under a uniform distribution
    over microstates in each macrocell:
      M[m', m] = fraction of microstates in macro m that map to macro m'
    Returns (macro_labels, M) where columns sum to 1.
    """
    states = enumerate_states(N)
    mac_labels = [coarse(s) for s in states]
    unique_macros = sorted(set(mac_labels), key=lambda x: (str(type(x)), str(x)))
    idx_of_macro = {m: i for i, m in enumerate(unique_macros)}
    # Collect microstates by macrocell
    by_macro: Dict[object, List[State]] = defaultdict(list)
    for s, m in zip(states, mac_labels):
        by_macro[m].append(s)
    # Build column-stochastic matrix
    M = np.zeros((len(unique_macros), len(unique_macros)), dtype=float)
    for m in unique_macros:
        bucket = by_macro[m]
        if not bucket:
            continue
        out_counts = Counter()
        for s in bucket:
            s_next = step_fn(s)
            m_next = coarse(s_next)
            out_counts[m_next] += 1
        col = np.zeros(len(unique_macros), dtype=float)
        denom = float(len(bucket))
        for m_next, c in out_counts.items():
            col[idx_of_macro[m_next]] = c / denom
        M[:, idx_of_macro[m]] = col
    return unique_macros, M


def is_markovianly_closed(
    N: int,
    coarse: CoarseFn,
    step_fn: StepFn,
) -> bool:
    """
    Test lumpability: for each macro m, all microstates in m induce the same
    distribution over next macrostates. Since dynamics is deterministic, this
    is equivalent to: for any s1, s2 in same macrocell, coarse(T(s1)) == coarse(T(s2)).
    """
    buckets: Dict[object, List[State]] = defaultdict(list)
    for s in enumerate_states(N):
        buckets[coarse(s)].append(s)
    for m, bucket in buckets.items():
        if len(bucket) <= 1:
            continue
        image = {str(coarse(step_fn(s))) for s in bucket}
        if len(image) != 1:
            return False
    return True


# ----------------------------
# Convenience experiments
# ----------------------------


def run_demo(
    N: int,
    steps: int,
    cg_func: CoarseFn,
    rule_str: str,
    eca_rule_num: int | None = None,
    layout: str = "spring",
) -> None:
    """Small demo printing macro entropy over time and visualizing trajectory on coarse-grained graph.

    Args:
        N: Number of bits
        steps: Number of evolution steps
        cg_func: Custom coarse-graining function (required)
        rule_str: Space-separated integers 0…2^N−1 representing the mapping (required)
        eca_rule_num: If provided, display 'Rule-<num>' instead of the full rule string
    """
    cg = cg_func
    step_fn = parse_permutation_rule(rule_str, N)

    s0 = int_to_state(1, N)  # start from 00..01
    traj = evolve(s0, steps, step_fn)
    macs = [cg(s) for s in traj]

    # Build macro-to-micros mapping to get entropy of each macrostate
    states = enumerate_states(N)
    macro_to_micros: Dict[object, List[int]] = defaultdict(list)
    for s in states:
        macro = cg(s)
        macro_to_micros[macro].append(state_to_int(s))

    # Compute entropy at each time step (entropy of the macrostate occupied)
    # This represents the uncertainty about which microstate we're in, given the macrostate
    entropies = []
    for mac in macs:
        num_micros = len(macro_to_micros[mac])
        # Entropy of uniform distribution over microstates in this macrostate
        H_t = math.log(num_micros) if num_micros > 0 else 0.0
        entropies.append(H_t)

    counts = Counter(macs)
    H_macro = entropy_from_counts(counts)  # single bag over the whole run

    rule_label = (
        f"Rule-{eca_rule_num}" if eca_rule_num is not None else f"Rule = {rule_str}"
    )
    print(f"N={N}, steps={steps}, {rule_label}")
    print(f"First 10 macrostates along trajectory: {macs[:10]}")
    print(f"Histogram of macrostates on trajectory: {dict(counts)}")
    print(f"Macro entropy over trajectory bag (nats): {H_macro:.6f}")
    print(f"Markovianly closed? {is_markovianly_closed(N, cg, step_fn)}")

    labels, M = induced_macro_transition_uniform(N, cg, step_fn)
    print("Induced macro transition matrix M[m', m] with uniform-in-macro weighting:")
    for i, row in enumerate(M):
        print(f"to {labels[i]}: " + " ".join(f"{x:.3f}" for x in row))

    # Visualize the trajectory on the coarse-grained graph
    _visualize_demo_trajectory(
        N, cg, macs, entropies, step_fn, rule_str, eca_rule_num, layout=layout
    )


def _visualize_demo_trajectory(
    N: int,
    cg: CoarseFn,
    trajectory_macs: List[object],
    entropies: List[float],
    step_fn: StepFn,
    rule: str,
    eca_rule_num: int | None = None,
    layout: str = "spring",
    coarse_node_size: int = 1200,
) -> None:
    """
    Visualize the trajectory on the coarse-grained phase space graph with entropy plot.
    """
    # Get all microstates and their macrostates
    states = enumerate_states(N)
    micro_to_macro = {state_to_int(s): cg(s) for s in states}

    # Group microstates by macrostate
    macro_to_micros: Dict[object, List[int]] = defaultdict(list)
    for micro_int, macro in micro_to_macro.items():
        macro_to_micros[macro].append(micro_int)

    # Get unique macrostates
    unique_macros = sorted(macro_to_micros.keys(), key=lambda x: (str(type(x)), str(x)))

    # Build macro transition graph
    G = nx.DiGraph()

    # Add nodes
    for macro in unique_macros:
        count = len(macro_to_micros[macro])
        if isinstance(macro, tuple):
            label = f"{''.join(map(str, macro))}\n({count})"
        else:
            label = f"{macro}\n({count})"
        G.add_node(macro, label=label, count=count)

    # Build macro transitions
    macro_transitions: Dict[Tuple[object, object], int] = defaultdict(int)
    for macro_from in unique_macros:
        for micro_int in macro_to_micros[macro_from]:
            s = int_to_state(micro_int, N)
            s_next = step_fn(s)
            macro_to = cg(s_next)
            macro_transitions[(macro_from, macro_to)] += 1

    # Add edges
    for (macro_from, macro_to), weight in macro_transitions.items():
        G.add_edge(macro_from, macro_to, weight=weight)

    # Create figure with two subplots: graph and entropy
    fig = plt.figure(figsize=(16, 8))
    gs = fig.add_gridspec(1, 2, width_ratios=[2, 1])
    ax_graph = fig.add_subplot(gs[0])
    ax_entropy = fig.add_subplot(gs[1])

    # Layout
    pos = _compute_layout(G, layout)

    # Node sizes: enforce a uniform size independent of microstate counts
    # This responds to the requirement that node size in the coarse-grained description
    # should NOT reflect the number of microstates. We keep counts only in labels.
    uniform_size = int(coarse_node_size)
    node_sizes = [uniform_size for _ in G.nodes()]

    # Draw base graph
    nx.draw_networkx_nodes(
        G, pos, node_color="lightcoral", node_size=node_sizes, alpha=0.3, ax=ax_graph
    )

    # Draw edges
    nx.draw_networkx_edges(
        G,
        pos,
        edge_color="gray",
        arrows=True,
        arrowsize=20,
        arrowstyle="->",
        width=2.5,
        ax=ax_graph,
        connectionstyle="arc3,rad=0.1",
        alpha=0.3,
    )

    # Draw labels
    labels = {macro: G.nodes[macro]["label"] for macro in G.nodes()}
    nx.draw_networkx_labels(
        G, pos, labels, font_size=9, font_weight="bold", ax=ax_graph
    )

    # Highlight the trajectory path
    traj_edges = []
    for i in range(len(trajectory_macs) - 1):
        m_from = trajectory_macs[i]
        m_to = trajectory_macs[i + 1]
        if G.has_edge(m_from, m_to):
            traj_edges.append((m_from, m_to))

    if traj_edges:
        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=traj_edges,
            edge_color="red",
            arrows=True,
            arrowsize=25,
            arrowstyle="->",
            width=3,
            ax=ax_graph,
            connectionstyle="arc3,rad=0.1",
        )

    # Highlight visited nodes
    visited_macros = list(set(trajectory_macs))
    # Visited nodes use the same uniform size
    visited_node_sizes = [uniform_size for _ in visited_macros]

    visited_pos = {m: pos[m] for m in visited_macros if m in pos}
    nx.draw_networkx_nodes(
        G,
        visited_pos,
        nodelist=visited_macros,
        node_color="darkred",
        node_size=visited_node_sizes,
        alpha=0.8,
        ax=ax_graph,
    )

    # Highlight starting position
    if trajectory_macs[0] in pos:
        start_macro = trajectory_macs[0]
        nx.draw_networkx_nodes(
            G,
            {start_macro: pos[start_macro]},
            nodelist=[start_macro],
            node_color="green",
            node_size=[uniform_size],
            alpha=1.0,
            ax=ax_graph,
            edgecolors="darkgreen",
            linewidths=3,
        )

    ax_graph.set_title(
        f"Trajectory Visualization (N={N}, {f'Rule-{eca_rule_num}' if eca_rule_num is not None else f'Rule = {rule}'})\n"
        f"Green=Start, Red edges=Trajectory path",
        fontsize=12,
        fontweight="bold",
    )
    ax_graph.axis("off")

    # Plot entropy over time
    time_steps = list(range(len(entropies)))
    ax_entropy.plot(time_steps, entropies, "b-o", linewidth=2, markersize=4)
    ax_entropy.set_xlabel("Time Step", fontsize=11, fontweight="bold")
    ax_entropy.set_ylabel("Macrostate Entropy (nats)", fontsize=11, fontweight="bold")
    ax_entropy.set_title("Macrostate Entropy Over Time", fontsize=12, fontweight="bold")
    ax_entropy.grid(True, alpha=0.3)
    ax_entropy.set_xlim(-0.5, len(entropies) - 0.5)

    plt.tight_layout()
    plt.show()


def print_cycles(N: int, rule_str: str, eca_rule_num: int | None = None) -> None:
    step_fn = parse_permutation_rule(rule_str, N)
    cycs = cycles(N, step_fn)
    Ls = sorted(len(c) for c in cycs)
    rule_label = (
        f"Rule-{eca_rule_num}" if eca_rule_num is not None else f"Rule = {rule_str}"
    )
    print(
        f"{rule_label}: Permutation decomposes into {len(cycs)} cycles. Lengths: {Ls}"
    )
    if N <= 5:
        for c in cycs:
            print([bin(x)[2:].zfill(N) for x in c])


def visualize_graph(
    N: int, rule_str: str, eca_rule_num: int | None = None, layout: str = "spring"
) -> None:
    """
    Visualize the directed graph of the microscopic phase space.
    Nodes are bit strings and edges connect each configuration to its successor under T.
    If eca_rule_num is provided, display 'Rule-<num>' instead of the full rule string.
    """
    step_fn = parse_permutation_rule(rule_str, N)
    mapping = permutation_graph(N, step_fn)

    # Create directed graph
    G = nx.DiGraph()

    # Add nodes with bit string labels
    for x in range(2**N):
        bit_string = bin(x)[2:].zfill(N)
        G.add_node(x, label=bit_string)

    # Add edges
    for x, y in mapping.items():
        G.add_edge(x, y)

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 10))

    pos = _compute_layout(G, layout)

    # Draw the graph
    nx.draw_networkx_nodes(G, pos, node_color="lightblue", node_size=800, ax=ax)
    nx.draw_networkx_edges(
        G,
        pos,
        edge_color="gray",
        arrows=True,
        arrowsize=20,
        arrowstyle="->",
        ax=ax,
        connectionstyle="arc3,rad=0.1",
    )

    # Draw labels as bit strings
    labels = {x: bin(x)[2:].zfill(N) for x in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, font_size=10, font_weight="bold", ax=ax)

    # Title: show Rule-<num> if eca_rule_num is provided
    if eca_rule_num is not None:
        rule_label = f"Rule-{eca_rule_num}"
    else:
        rule_label = "Rule = " + rule_str
    ax.set_title(
        f"Phase Portrait Graph (N={N}, {rule_label})\n$|\\Omega| = 2^{{{N}}} = {2**N}$ states",
        fontsize=14,
        fontweight="bold",
    )
    ax.axis("off")
    plt.tight_layout()

    plt.show()

    # Print cycle information
    cycs = cycles(N, step_fn)
    print("\nGraph structure:")
    print(f"  Nodes: {G.number_of_nodes()}")
    print(f"  Edges: {G.number_of_edges()}")
    print(f"  Cycles: {len(cycs)}")
    print(f"  Cycle lengths: {sorted(len(c) for c in cycs)}")


def visualize_coarse_graph(
    N: int,
    cg_func: CoarseFn,
    rule_str: str,
    eca_rule_num: int | None = None,
    layout: str = "spring",
    coarse_node_size: int = 1200,
) -> None:
    """
    Visualize the coarse-grained phase space graph.
    Nodes are macrostates (labeled with their macrostate value and microstate count),
    and edges show transitions between macrostates.

    Args:
        N: Number of bits
        cg_func: Custom coarse-graining function (required)
        rule_str: Space-separated integers representing the mapping
        eca_rule_num: If provided, display 'Rule-<num>' instead of the full rule string
    """
    if cg_func is None:
        raise ValueError(
            "Coarse-graining function is required. Use --groups to specify a custom partition."
        )
    cg = cg_func
    step_fn = parse_permutation_rule(rule_str, N)

    # Get all microstates and their macrostates
    states = enumerate_states(N)
    micro_to_macro = {state_to_int(s): cg(s) for s in states}

    # Group microstates by macrostate
    macro_to_micros: Dict[object, List[int]] = defaultdict(list)
    for micro_int, macro in micro_to_macro.items():
        macro_to_micros[macro].append(micro_int)

    # Get unique macrostates (sorted for consistent ordering)
    unique_macros = sorted(macro_to_micros.keys(), key=lambda x: (str(type(x)), str(x)))

    # Build macro transition graph
    G = nx.DiGraph()

    # Add nodes for each macrostate with label showing macro value and microstate count
    for macro in unique_macros:
        count = len(macro_to_micros[macro])
        # Format the label based on the type of macrostate
        if isinstance(macro, tuple):
            # For rotation class, show the bit pattern
            label = f"{''.join(map(str, macro))}\n({count})"
        else:
            label = f"{macro}\n({count})"
        G.add_node(macro, label=label, count=count)

    # Build macro transitions: for each macrostate, see where its microstates go
    macro_transitions: Dict[Tuple[object, object], int] = defaultdict(int)
    for macro_from in unique_macros:
        for micro_int in macro_to_micros[macro_from]:
            s = int_to_state(micro_int, N)
            s_next = step_fn(s)
            macro_to = cg(s_next)
            macro_transitions[(macro_from, macro_to)] += 1

    # Add edges with weights (number of microstates making this transition)
    for (macro_from, macro_to), weight in macro_transitions.items():
        G.add_edge(macro_from, macro_to, weight=weight)

    # Create figure
    fig, ax = plt.subplots(figsize=(14, 10))

    pos = _compute_layout(G, layout)

    uniform_size = int(coarse_node_size)
    node_sizes = [uniform_size for _ in G.nodes()]

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos, node_color="lightcoral", node_size=node_sizes, alpha=0.7, ax=ax
    )

    nx.draw_networkx_edges(
        G,
        pos,
        edge_color="gray",
        arrows=True,
        arrowsize=20,
        arrowstyle="->",
        width=2.5,
        ax=ax,
        connectionstyle="arc3,rad=0.1",
    )

    # Draw labels
    labels = {macro: G.nodes[macro]["label"] for macro in G.nodes()}
    nx.draw_networkx_labels(G, pos, labels, font_size=9, font_weight="bold", ax=ax)

    # Always display edge labels showing transition counts
    edge_labels = {
        (u, v): f"{G[u][v]['weight']}" for u, v in G.edges() if G[u][v]["weight"] > 0
    }
    nx.draw_networkx_edge_labels(
        G,
        pos,
        edge_labels,
        font_size=10,
        font_color="blue",
        ax=ax,
        connectionstyle="arc3,rad=0.1",
    )

    # Title with coarse-graining info
    rule_label = (
        f"Rule-{eca_rule_num}" if eca_rule_num is not None else f"Rule = {rule_str}"
    )
    ax.set_title(
        f"Coarse-Grained Phase Portrait Graph (N={N}, {rule_label})\n"
        f"Macrostates: {len(unique_macros)}, Total microstates: {2**N}\n",
        fontsize=12,
        fontweight="bold",
    )
    ax.axis("off")
    plt.tight_layout()

    plt.show()

    # Print statistics
    print("\nCoarse-grained graph structure:")
    print(f"  Macrostates: {len(unique_macros)}")
    print("  Microstates per macrostate:")
    for macro in unique_macros:
        count = len(macro_to_micros[macro])
        print(f"    {macro}: {count} microstates")
    print(f"  Macro transitions: {G.number_of_edges()}")
    print(f"  Markovianly closed? {is_markovianly_closed(N, cg, step_fn)}")


# ----------------------------
# CLI
# ----------------------------


def _parse_custom_partition(spec: str, N: int) -> CoarseFn:
    """Parse a user supplied partition of the 2**N microstates and return a coarse function.

    The *spec* can be either:
      1. A path to a JSON file
      2. An inline JSON string

    Accepted JSON formats:
      - Mapping form: {"labelA": ["000", 3, "5"], "labelB": [...], ...}
      - List form: [["000", 1], ["010", "011"], ...]  (labels become "0", "1", ...)

    Microstates inside groups can be:
      * Bit strings of length N (e.g. "0101")
      * Integers (e.g. 5) or numeric strings ("5") representing the integer encoding

    Validation ensures the groups form a *partition*: every of the 2**N states
    appears in exactly one group. Duplicate or missing states raise ValueError.
    """
    # Load text
    if os.path.exists(spec):
        with open(spec, "r", encoding="utf-8") as f:
            text = f.read()
    else:
        text = spec
    try:
        data = json.loads(text)
    except json.JSONDecodeError as e:
        raise ValueError(
            f"Failed to parse custom grouping JSON (position {e.pos}): {e.msg}"
        ) from e

    # Normalize to list of (label, list_of_state_ints)
    groups: List[Tuple[str, List[int]]] = []
    if isinstance(data, dict):
        for label, members in data.items():
            if not isinstance(members, list):
                raise ValueError(
                    f"Group '{label}' must map to a list, got {type(members)}"
                )
            groups.append((str(label), members))
    elif isinstance(data, list):
        for i, members in enumerate(data):
            if not isinstance(members, list):
                raise ValueError(
                    f"Element {i} of top-level list must be a list (group), got {type(members)}"
                )
            groups.append((str(i), members))
    else:
        raise ValueError(
            "Custom grouping JSON must be dict or list. See README for accepted formats."
        )

    total_states = 2**N
    all_states_set = set(range(total_states))
    seen: Dict[int, str] = {}
    int_to_label: Dict[int, str] = {}

    def parse_member(m) -> int:
        if isinstance(m, int):
            return m
        if isinstance(m, str):
            m_str = m.strip()
            if set(m_str).issubset({"0", "1"}) and len(m_str) == N:
                # Bitstring form
                val = 0
                for ch in m_str:
                    val = (val << 1) | int(ch)
                return val
            # Try as integer literal
            try:
                return int(m_str)
            except ValueError:
                raise ValueError(
                    f"Cannot interpret member '{m}' as bitstring of length {N} or integer"
                )
        raise ValueError(f"Unsupported member type: {type(m)} for {m}")

    for label, members in groups:
        for m in members:
            val = parse_member(m)
            if not (0 <= val < total_states):
                raise ValueError(
                    f"State value {val} (from {m}) out of range for N={N} (0..{total_states - 1})"
                )
            if val in seen:
                raise ValueError(
                    f"Duplicate assignment: state {val} already in group '{seen[val]}', cannot also be in '{label}'"
                )
            seen[val] = label
            int_to_label[val] = label

    missing = all_states_set.difference(seen.keys())
    if missing:
        # Provide a helpful hint by showing a few missing states in both int and bitstring forms
        sample = sorted(list(missing))[:5]

        def fmt(v: int) -> str:
            return bin(v)[2:].zfill(N)

        raise ValueError(
            "Custom grouping is not a full partition. Missing states: "
            + ", ".join(f"{v} ('{fmt(v)}')" for v in sample)
            + ("..." if len(missing) > 5 else "")
        )

    # Build coarse function closure
    def custom_cg(s: State) -> str:
        idx = state_to_int(s)
        return int_to_label[idx]

    return custom_cg


def main():
    parser = argparse.ArgumentParser(
        description="Finite Descriptive Universe simulator (general reversible rules)."
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_demo = sub.add_parser("demo", help="Run a small simulation and print macro info.")
    p_demo.add_argument("-N", type=int, default=4, help="Number of bits")
    p_demo.add_argument("-t", "--steps", type=int, default=16, help="Number of steps")
    demo_rule_group = p_demo.add_mutually_exclusive_group(required=True)
    demo_rule_group.add_argument(
        "-r",
        "--rule",
        type=str,
        help="Space-separated integers 0…2^N−1 representing the mapping from each microstate to its image (e.g. '2 3 0 1' for N=2)",
    )
    demo_rule_group.add_argument(
        "--eca",
        type=int,
        metavar="RULE_NUM",
        help="Elementary cellular automaton rule number (0-255). Automatically generates the full 2^N mapping.",
    )
    p_demo.add_argument(
        "--groups",
        type=str,
        required=True,
        help="Custom grouping specification (JSON string or path). Required.",
    )
    p_demo.add_argument(
        "--plot",
        action="store_true",
        help="Show coarse-grained phase space graph",
    )
    # Layout option (shared with graph visual commands)
    p_demo.add_argument(
        "--layout",
        type=str,
        default="spring",
        help="Layout algorithm for trajectory coarse graph when --plot is used (spring, circular, shell, kamada, spectral, planar).",
    )
    p_demo.add_argument(
        "--coarse-node-size",
        type=int,
        default=1200,
        help="Uniform node size for coarse-grained graph nodes (trajectory & coarse-graph). Counts appear only in labels.",
    )

    p_cycles = sub.add_parser("cycles", help="Print cycle decomposition stats.")
    p_cycles.add_argument("-N", type=int, default=4, help="Number of bits")
    cycles_rule_group = p_cycles.add_mutually_exclusive_group(required=True)
    cycles_rule_group.add_argument(
        "-r",
        "--rule",
        type=str,
        help="Space-separated integers 0…2^N−1 representing the mapping from each microstate to its image",
    )
    cycles_rule_group.add_argument(
        "--eca",
        type=int,
        metavar="RULE_NUM",
        help="Elementary cellular automaton rule number (0-255). Automatically generates the full 2^N mapping.",
    )

    p_graph = sub.add_parser(
        "graph", help="Visualize the microscopic phase space graph."
    )
    p_graph.add_argument("-N", type=int, default=4, help="Number of bits")
    p_graph.add_argument(
        "--layout",
        type=str,
        default="spring",
        help="Graph layout algorithm (spring, circular, shell, kamada, spectral, planar).",
    )
    graph_rule_group = p_graph.add_mutually_exclusive_group(required=True)
    graph_rule_group.add_argument(
        "-r",
        "--rule",
        type=str,
        help="Space-separated integers 0…2^N−1 representing the mapping from each microstate to its image",
    )
    graph_rule_group.add_argument(
        "--eca",
        type=int,
        metavar="RULE_NUM",
        help="Elementary cellular automaton rule number (0-255). Automatically generates the full 2^N mapping.",
    )

    p_coarse_graph = sub.add_parser(
        "coarse-graph", help="Visualize the coarse-grained phase space graph."
    )
    p_coarse_graph.add_argument("-N", type=int, default=4, help="Number of bits")
    p_coarse_graph.add_argument(
        "--layout",
        type=str,
        default="spring",
        help="Graph layout algorithm (spring, circular, shell, kamada, spectral, planar).",
    )
    coarse_rule_group = p_coarse_graph.add_mutually_exclusive_group(required=True)
    coarse_rule_group.add_argument(
        "-r",
        "--rule",
        type=str,
        help="Space-separated integers 0…2^N−1 representing the mapping from each microstate to its image",
    )
    coarse_rule_group.add_argument(
        "--eca",
        type=int,
        metavar="RULE_NUM",
        help="Elementary cellular automaton rule number (0-255). Automatically generates the full 2^N mapping.",
    )
    p_coarse_graph.add_argument(
        "--groups",
        type=str,
        required=True,
        help="Custom grouping specification (JSON string or path). Required.",
    )
    p_coarse_graph.add_argument(
        "--coarse-node-size",
        type=int,
        default=1200,
        help="Uniform node size for coarse-grained graph nodes. Overrides prior proportional sizing.",
    )

    args = parser.parse_args()

    # Convert --eca to rule string if provided
    if hasattr(args, "eca") and args.eca is not None:
        args.rule = generate_eca_rule_string(args.eca, args.N)

    if args.cmd == "demo":
        cg_func = _parse_custom_partition(args.groups, args.N)
        run_demo(args.N, args.steps, cg_func, args.rule, args.eca, layout=args.layout)
        if args.plot:
            visualize_coarse_graph(
                args.N,
                cg_func,
                args.rule,
                args.eca,
                layout=args.layout,
                coarse_node_size=args.coarse_node_size,
            )
    elif args.cmd == "cycles":
        print_cycles(args.N, args.rule, args.eca)
    elif args.cmd == "graph":
        visualize_graph(args.N, args.rule, args.eca, layout=args.layout)
    elif args.cmd == "coarse-graph":
        cg_func = _parse_custom_partition(args.groups, args.N)
        visualize_coarse_graph(
            args.N,
            cg_func,
            args.rule,
            args.eca,
            layout=args.layout,
            coarse_node_size=args.coarse_node_size,
        )


def _compute_layout(G: nx.Graph, layout: str) -> Dict[object, np.ndarray]:
    """Compute node positions for a requested layout algorithm.

    Supported layout strings (case-insensitive): spring, circular, shell, kamada,
    spectral, planar. Falls back to spring if the requested layout is invalid
    or planar layout fails (e.g., graph not planar).
    """
    key = (layout or "spring").lower()
    try:
        if key == "spring":
            return nx.spring_layout(G, k=3, iterations=200, seed=42)
        if key == "circular":
            return nx.circular_layout(G)
        if key == "shell":
            return nx.shell_layout(G)
        if key == "kamada":
            return nx.kamada_kawai_layout(G)
        if key == "spectral":
            return nx.spectral_layout(G)
        if key == "planar":
            try:
                return nx.planar_layout(G)
            except nx.NetworkXException:
                pass
    except Exception:
        pass
    return nx.spring_layout(G, k=None, iterations=300, seed=42)


if __name__ == "__main__":
    main()
