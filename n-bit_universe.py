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
    """Small demo printing macro entropy over time and closure result."""
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
    counts = Counter(macs)
    H_macro = entropy_from_counts(counts)  # single bag over the whole run
    # If you want H_macro(t), compute a sliding histogram or per-time distribution.

    print(f"N={N}, steps={steps}, coarse={coarse}")
    print(f"First 10 macrostates along trajectory: {macs[:10]}")
    print(f"Histogram of macrostates on trajectory: {dict(counts)}")
    print(f"Macro entropy over trajectory bag (nats): {H_macro:.6f}")
    print(f"Markovianly closed? {is_markovianly_closed(N, cg, rule90_step)}")

    labels, M = induced_macro_transition_uniform(N, cg, rule90_step)
    print("Induced macro transition matrix M[m', m] with uniform-in-macro weighting:")
    for i, row in enumerate(M):
        print(f"to {labels[i]}: " + " ".join(f"{x:.3f}" for x in row))


def print_cycles(N: int) -> None:
    cycs = cycles(N, rule90_step)
    Ls = sorted(len(c) for c in cycs)
    print(f"Permutation decomposes into {len(cycs)} cycles. Lengths: {Ls}")
    # Optionally print the cycles themselves for very small N
    if N <= 5:
        for c in cycs:
            print([bin(x)[2:].zfill(N) for x in c])


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

    p_cycles = sub.add_parser("cycles", help="Print cycle decomposition stats.")
    p_cycles.add_argument("-N", type=int, default=4, help="Number of bits")

    args = parser.parse_args()
    if args.cmd == "demo":
        run_demo(args.N, args.steps, args.coarse)
    elif args.cmd == "cycles":
        print_cycles(args.N)


if __name__ == "__main__":
    main()
