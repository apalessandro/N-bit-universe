"""
Finite Descriptive Universe: rule-90 cellular automaton on N bits.
Implements:
  - Deterministic synchronous updates with periodic boundary conditions
  - Coarse-grainings: Hamming weight, parity, rotation class
  - Microscopic and macroscopic entropy (natural logs)
  - Induced macro transition operator (uniform over microstates within a macrocell)
  - Test for Markovian closure (lumpability)
  - Cycle decomposition of the microscopic permutation
  - Simple CLI for quick experiments
"""

from __future__ import annotations
import argparse
import math
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


def rule90_step(s: State) -> State:
    """One synchronous update of rule-90 with periodic boundaries."""
    return np.bitwise_xor(np.roll(s, 1), np.roll(s, -1)).astype(np.uint8)


def evolve(s0: State, steps: int, step_fn: StepFn = rule90_step) -> List[State]:
    """Return the trajectory [s0, s1, ..., s_steps]."""
    traj = [s0.copy()]
    s = s0.copy()
    for _ in range(steps):
        s = step_fn(s)
        traj.append(s.copy())
    return traj


# ----------------------------
# Coarse-grainings
# ----------------------------


def cg_weight(s: State) -> int:
    """Hamming weight coarse-graining."""
    return int(s.sum())


def cg_parity(s: State) -> int:
    """Parity coarse-graining, 0 even, 1 odd."""
    return int(s.sum() & 1)


def _min_rotation(bits: Tuple[int, ...]) -> Tuple[int, ...]:
    """Minimal lexicographic rotation of the bit tuple."""
    N = len(bits)
    doubles = bits + bits
    # Booth algorithm not needed for small N
    rotations = (tuple(doubles[i : i + N]) for i in range(N))
    return min(rotations)


def cg_rotation_class(s: State) -> Tuple[int, ...]:
    """Rotation class representative as the minimal lexicographic rotation."""
    return _min_rotation(tuple(int(x) for x in s.tolist()))


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


def permutation_graph(N: int, step_fn: StepFn = rule90_step) -> Dict[int, int]:
    """Return the permutation as a dict: x -> y in integer labels."""
    mapping: Dict[int, int] = {}
    for s in enumerate_states(N):
        x = state_to_int(s)
        y = state_to_int(step_fn(s))
        mapping[x] = y
    return mapping


def cycles(N: int, step_fn: StepFn = rule90_step) -> List[List[int]]:
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
    step_fn: StepFn = rule90_step,
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
    step_fn: StepFn = rule90_step,
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


def run_demo(N: int, steps: int, coarse: str) -> None:
    """Small demo printing macro entropy over time and visualizing trajectory on coarse-grained graph."""
    coarse_map: Dict[str, CoarseFn] = {
        "weight": cg_weight,
        "parity": cg_parity,
        "rotation": cg_rotation_class,
    }
    if coarse not in coarse_map:
        raise ValueError(
            f"Unknown coarse-graining {coarse}. Choose from {list(coarse_map)}"
        )
    cg = coarse_map[coarse]

    s0 = int_to_state(1, N)  # start from 00..01
    traj = evolve(s0, steps, rule90_step)
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

    print(f"N={N}, steps={steps}, coarse={coarse}")
    print(f"First 10 macrostates along trajectory: {macs[:10]}")
    print(f"Histogram of macrostates on trajectory: {dict(counts)}")
    print(f"Macro entropy over trajectory bag (nats): {H_macro:.6f}")
    print(f"Markovianly closed? {is_markovianly_closed(N, cg, rule90_step)}")

    labels, M = induced_macro_transition_uniform(N, cg, rule90_step)
    print("Induced macro transition matrix M[m', m] with uniform-in-macro weighting:")
    for i, row in enumerate(M):
        print(f"to {labels[i]}: " + " ".join(f"{x:.3f}" for x in row))

    # Visualize the trajectory on the coarse-grained graph
    _visualize_demo_trajectory(N, cg, macs, entropies, coarse, rule90_step)


def _visualize_demo_trajectory(
    N: int,
    cg: CoarseFn,
    trajectory_macs: List[object],
    entropies: List[float],
    coarse: str,
    step_fn: StepFn = rule90_step,
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
    if len(unique_macros) <= 10:
        pos = nx.spring_layout(G, k=3, iterations=50, seed=42)
    else:
        pos = nx.kamada_kawai_layout(G)

    # Node sizes
    max_node_size = 3000
    raw_sizes = [G.nodes[macro]["count"] * 300 for macro in G.nodes()]
    node_sizes = [min(size, max_node_size) for size in raw_sizes]

    # Draw base graph
    nx.draw_networkx_nodes(
        G, pos, node_color="lightcoral", node_size=node_sizes, alpha=0.3, ax=ax_graph
    )

    # Draw edges
    edges = G.edges()
    weights = [G[u][v]["weight"] for u, v in edges]
    max_weight = max(weights) if weights else 1
    edge_widths = [2 + 3 * (w / max_weight) for w in weights]

    nx.draw_networkx_edges(
        G,
        pos,
        edge_color="gray",
        arrows=True,
        arrowsize=20,
        arrowstyle="->",
        width=edge_widths,
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
    visited_node_sizes = []
    for macro in visited_macros:
        if macro in G.nodes():
            idx = list(G.nodes()).index(macro)
            visited_node_sizes.append(node_sizes[idx])
        else:
            visited_node_sizes.append(800)

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
            node_size=[node_sizes[list(G.nodes()).index(start_macro)]],
            alpha=1.0,
            ax=ax_graph,
            edgecolors="darkgreen",
            linewidths=3,
        )

    coarse_display = coarse.capitalize()
    ax_graph.set_title(
        f"Trajectory Visualization (N={N}, {coarse_display})\n"
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


def print_cycles(N: int) -> None:
    cycs = cycles(N, rule90_step)
    Ls = sorted(len(c) for c in cycs)
    print(f"Permutation decomposes into {len(cycs)} cycles. Lengths: {Ls}")
    # Optionally print the cycles themselves for very small N
    if N <= 5:
        for c in cycs:
            print([bin(x)[2:].zfill(N) for x in c])


def visualize_graph(N: int, step_fn: StepFn = rule90_step) -> None:
    """
    Visualize the directed graph of the microscopic phase space.
    Nodes are bit strings and edges connect each configuration to its successor under T.
    """
    # Get the permutation graph
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

    # Use circular layout for better cycle visualization
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

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

    ax.set_title(
        f"Microscopic Phase Space Graph for N={N}\n"
        f"$|\\Omega| = 2^{{{N}}} = {2**N}$ states",
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
    coarse: str,
    step_fn: StepFn = rule90_step,
) -> None:
    """
    Visualize the coarse-grained phase space graph.
    Nodes are macrostates (labeled with their macrostate value and microstate count),
    and edges show transitions between macrostates.
    """
    coarse_map: Dict[str, CoarseFn] = {
        "weight": cg_weight,
        "parity": cg_parity,
        "rotation": cg_rotation_class,
    }
    if coarse not in coarse_map:
        raise ValueError(
            f"Unknown coarse-graining {coarse}. Choose from {list(coarse_map)}"
        )
    cg = coarse_map[coarse]

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

    # Use layout appropriate for the graph structure
    if len(unique_macros) <= 10:
        pos = nx.spring_layout(G, k=3, iterations=50, seed=42)
    else:
        pos = nx.kamada_kawai_layout(G)

    # Node sizes proportional to number of microstates with a maximum cap
    max_node_size = 3000  # Maximum size for the red ball
    raw_sizes = [G.nodes[macro]["count"] * 300 for macro in G.nodes()]
    node_sizes = [min(size, max_node_size) for size in raw_sizes]

    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos, node_color="lightcoral", node_size=node_sizes, alpha=0.7, ax=ax
    )

    # Draw edges with varying width based on weight
    edges = G.edges()
    weights = [G[u][v]["weight"] for u, v in edges]
    max_weight = max(weights) if weights else 1
    edge_widths = [2 + 3 * (w / max_weight) for w in weights]

    nx.draw_networkx_edges(
        G,
        pos,
        edge_color="gray",
        arrows=True,
        arrowsize=20,
        arrowstyle="->",
        width=edge_widths,
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
        G, pos, edge_labels, font_size=10, font_color="blue", ax=ax
    )

    # Title with coarse-graining info
    coarse_display = coarse.capitalize()
    ax.set_title(
        f"Coarse-Grained Phase Space Graph (N={N}, Coarse-graining: {coarse_display})\n"
        f"Macrostates: {len(unique_macros)}, Total microstates: {2**N}\n"
        f"Node labels show: macrostate value (microstate count)",
        fontsize=12,
        fontweight="bold",
    )
    ax.axis("off")
    plt.tight_layout()

    plt.show()

    # Print statistics
    print("\nCoarse-grained graph structure:")
    print(f"  Coarse-graining: {coarse_display}")
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


def main():
    parser = argparse.ArgumentParser(
        description="Finite Descriptive Universe simulator (rule-90)."
    )
    sub = parser.add_subparsers(dest="cmd", required=True)

    p_demo = sub.add_parser("demo", help="Run a small simulation and print macro info.")
    p_demo.add_argument("-N", type=int, default=4, help="Number of bits")
    p_demo.add_argument("-t", "--steps", type=int, default=16, help="Number of steps")
    p_demo.add_argument(
        "-c",
        "--coarse",
        type=str,
        default="parity",
        choices=["parity", "weight", "rotation"],
        help="Coarse-graining",
    )
    p_demo.add_argument(
        "--plot",
        action="store_true",
        help="Show coarse-grained phase space graph",
    )

    p_cycles = sub.add_parser("cycles", help="Print cycle decomposition stats.")
    p_cycles.add_argument("-N", type=int, default=4, help="Number of bits")

    p_graph = sub.add_parser(
        "graph", help="Visualize the microscopic phase space graph."
    )
    p_graph.add_argument("-N", type=int, default=4, help="Number of bits")

    p_coarse_graph = sub.add_parser(
        "coarse-graph", help="Visualize the coarse-grained phase space graph."
    )
    p_coarse_graph.add_argument("-N", type=int, default=4, help="Number of bits")
    p_coarse_graph.add_argument(
        "-c",
        "--coarse",
        type=str,
        default="parity",
        choices=["parity", "weight", "rotation"],
        help="Coarse-graining",
    )

    args = parser.parse_args()
    if args.cmd == "demo":
        run_demo(args.N, args.steps, args.coarse)
        if args.plot:
            visualize_coarse_graph(args.N, args.coarse)
    elif args.cmd == "cycles":
        print_cycles(args.N)
    elif args.cmd == "graph":
        visualize_graph(args.N)
    elif args.cmd == "coarse-graph":
        visualize_coarse_graph(args.N, args.coarse)


if __name__ == "__main__":
    main()
