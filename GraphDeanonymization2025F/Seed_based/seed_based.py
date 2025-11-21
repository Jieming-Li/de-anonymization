import argparse
from collections import defaultdict
import math

#read edgelist files and store graph in a dictionary
def read_edges(path):
    adj = defaultdict(set)
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            u, v = line.split()
            adj[u].add(v)
            adj[v].add(u)
    return adj

#read seeds and store in a dictionary
def read_seeds(path):
    seeds = {}
    with open(path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            u, v = line.split()
            seeds[u] = v
    return seeds

#write the final mapping
def write_mapping(mapping, path):
    with open(path, "w") as f:
        for u, v in mapping.items():
            f.write(f"{u} {v}\n")


#de-anonymization core

#compute score as: count of matched neighbors of u,v/sqrt(degree(u))*sqrt(degree(v))
def score_uv(u, v, adj1, adj2, matched_1_to_2):
    neighbors_u = adj1.get(u, ())
    neighbors_v = adj2.get(v, ())
    deg_u = len(neighbors_u)
    deg_v = len(neighbors_v)

    if deg_u == 0 or deg_v == 0:
        return 0.0

    matched_count = 0
    for nb in neighbors_u:
        if nb in matched_1_to_2:
            mapped_nb = matched_1_to_2[nb]
            if mapped_nb in neighbors_v:
                matched_count += 1

    denom = math.sqrt(deg_u) * math.sqrt(deg_v)
    return matched_count / denom if denom > 0 else 0.0

#iteratively expand seed mapping
def propagate(adj1, adj2, seeds, min_score, score_ratio, max_iters):
    matched_1_to_2 = dict(seeds)
    matched_2_to_1 = {v: u for u, v in seeds.items()}

    #get all nodes, including ones that are only neighbors
    nodes1 = set(adj1.keys())
    for s in adj1.values(): nodes1 |= set(s)

    nodes2 = set(adj2.keys())
    for s in adj2.values(): nodes2 |= set(s)

    unmatched1 = nodes1 - set(matched_1_to_2.keys())
    unmatched2 = nodes2 - set(matched_2_to_1.keys())

    for it in range(1, max_iters + 1):

        new_matches = {}

        for u in list(unmatched1):

            #build candidate set for u
            mapped_neighbors = []
            for nb in adj1.get(u, ()):
                if nb in matched_1_to_2:
                    mapped_neighbors.append(matched_1_to_2[nb])

            if not mapped_neighbors:
                continue

            #candidates: neighbors of the mapped neighbors
            candidates = set()
            for mv in mapped_neighbors:
                candidates.update(adj2.get(mv, ()))

            #remove already matched nodes
            candidates -= unmatched2.symmetric_difference(unmatched2 & candidates)

            if not candidates:
                continue

            #compute score for all candidates
            best = None
            best_score = -1
            second_score = -1

            for v in candidates:
                s = score_uv(u, v, adj1, adj2, matched_1_to_2)

                if s > best_score:
                    second_score = best_score
                    best = (v, s)
                    best_score = s
                elif s > second_score:
                    second_score = s

            if best is None:
                continue

            v_best, best_val = best
            sec_val = second_score

            accept = False

            #absolute threshold
            if best_val >= min_score:
                accept = True

            #ratio rule
            if sec_val == 0:
                accept = True
            else:
                if best_val / (sec_val + 1e-12) >= score_ratio:
                    accept = True

            if accept:
                new_matches[u] = v_best

        if not new_matches:
            break

        #commit new matches
        for u, v in new_matches.items():
            if v in matched_2_to_1:
                continue
            matched_1_to_2[u] = v
            matched_2_to_1[v] = u
            unmatched1.discard(u)
            unmatched2.discard(v)


    return matched_1_to_2


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--min_score", type=float, default=0.05,
                        help="Minimum normalized score to accept")
    parser.add_argument("--ratio", type=float, default=1.2,
                        help="Minimum ratio between best and second-best")
    parser.add_argument("--max_iters", type=int, default=100)
    args = parser.parse_args()

    adj1 = "seed_G1.edgelist"
    adj2 = "seed_G2.edgelist"
    seeds = "seed_mapping.txt"

    final_map = propagate(
        adj1, adj2, seeds,
        min_score=args.min_score,
        score_ratio=args.ratio,
        max_iters=args.max_iters
    )

    out="out.txt"
    write_mapping(final_map, out)

if __name__ == "__main__":
    main()
