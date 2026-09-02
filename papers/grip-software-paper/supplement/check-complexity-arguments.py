#!/usr/bin/env python3
"""Independent finite checks of S1 arguments, not a grip timing benchmark.

Uses exact integer-scaled path lengths and exact rational expectations.
Only the first genuine coarse level and its level-zero insertion queries
are modeled; this is not a replacement implementation of weighted GRIP.
"""

from fractions import Fraction
from functools import lru_cache
from heapq import heappop, heappush
from itertools import combinations
from math import isclose, log
from random import Random


def check_telescoping():
    count = 0
    for n in range(1, 13):
        for q in range(1, n + 1):
            interior = range(q + 1, n)
            for length in range(len(interior) + 1):
                for middle in combinations(interior, length):
                    levels = [n, *reversed(middle), q]
                    # Repeated levels contribute zero and must also be allowed.
                    levels.insert(1, n)
                    lhs = sum(1 - y / x for x, y in zip(levels, levels[1:]))
                    telescope = sum(log(x / y) for x, y in zip(levels, levels[1:]))
                    assert lhs <= log(n / q) + 1e-12
                    assert isclose(telescope, log(n / q), abs_tol=1e-12)
                    count += 1
    print(f"Telescoping inequality: {count} nested size sequences passed")


def check_repulsion():
    count = 0
    for population in range(3, 11):
        for sample in range(1, population - 1):
            @lru_cache(None)
            def remaining(accepted):
                if len(accepted) == sample:
                    return Fraction(0)
                new_partners = [v for v in range(1, population) if v not in accepted]
                # Solve the one-step recurrence, including self/duplicate draws.
                future = sum(
                    (remaining(accepted | {v}) for v in new_partners), Fraction(0)
                )
                return (population + future) / len(new_partners)

            formula = sum(
                (Fraction(population, population - 1 - k) for k in range(sample)),
                Fraction(0),
            )
            assert remaining(frozenset()) == formula
            assert formula <= Fraction(sample * population, population - sample)
            assert formula <= Fraction(sample * (sample + 2), 2)
            count += 1
    print(f"Repulsion expectation: {count} exact finite-state checks passed")


def traverse(adjacency, root, cutoff=None, stop=None):
    distances = {root: 0}
    heap = [(0, root)]
    prefix = []
    while heap:
        distance, vertex = heappop(heap)
        if distance != distances[vertex]:
            continue
        if cutoff is not None and distance > cutoff:
            break
        prefix.append(vertex)
        if vertex != root and stop is not None and stop(vertex):
            break  # The stopping vertex is settled but not expanded.
        for neighbor, weight in adjacency[vertex]:
            candidate = distance + weight
            if cutoff is not None and candidate > cutoff:
                continue
            if neighbor not in distances or candidate < distances[neighbor]:
                distances[neighbor] = candidate
                heappush(heap, (candidate, neighbor))
    return prefix


def check_compressed_path():
    print("Compressed path: n, order, coarse size, queries, min prefix, total prefix")
    for block_size in (96, 192, 384):
        n = 2 * block_size
        # Multiply normalized lengths by 4m: no floating-point boundary issues.
        scale = 4 * block_size
        lengths = [1] * (block_size - 1) + [scale] * block_size
        assert sorted(lengths)[len(lengths) // 2] == scale
        adjacency = [[] for _ in range(n)]
        for v, weight in enumerate(lengths):
            adjacency[v].append((v + 1, weight))
            adjacency[v + 1].append((v, weight))
        shuffled = list(range(n))
        Random(20260902).shuffle(shuffled)
        orders = {"forward": list(range(n)), "reverse": list(reversed(range(n))),
                  "shuffled": shuffled}
        block = set(range(block_size))
        for name, order in orders.items():
            candidates = set(range(n))
            coarse = set()
            for center in order:
                if center in candidates:
                    coarse.add(center)
                    candidates.difference_update(traverse(adjacency, center, cutoff=scale))
            assert not candidates
            assert len(coarse & block) <= 1
            assert len(coarse) >= block_size / 3
            assert len(coarse) > 24  # Genuine first level, not the top subset.
            sizes = []
            for root in sorted(block - coarse):
                cached = 0
                cache_complete = False
                anchors = 0

                def stop(vertex):
                    nonlocal cached, cache_complete, anchors
                    if cached < 12:
                        cached += 1
                    else:
                        cache_complete = True  # Extra eligible vertex, as in S1.
                    if vertex in coarse:
                        anchors += 1
                    return cache_complete and anchors >= 3

                prefix = traverse(adjacency, root, stop=stop)
                assert set(prefix[:block_size]) == block
                assert anchors >= 3 and cache_complete
                sizes.append(len(prefix))
            assert len(sizes) >= block_size - 1
            assert sum(sizes) >= block_size * (block_size - 1)
            print(n, name, len(coarse), len(sizes), min(sizes), sum(sizes), sep=", ")


if __name__ == "__main__":
    check_telescoping()
    check_repulsion()
    check_compressed_path()
    print("All checks passed; these are mathematical checks, not runtime measurements.")
