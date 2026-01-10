# Contig Classification Legend

This document explains the topological tags assigned to assembled contigs. Tags are generated based on the contig's relationship with the De Bruijn Graph structure (hubs, bubbles, tips, and cycles).

## Tagging Hierarchy
To ensure accuracy, tags are assigned following this priority:
1. **Tips** (Probable artifacts)
2. **Circular** (Closed loops)
3. **Lariat** (Lasso structures)
4. **Linear** (Standard paths)

---

## 1. Primary Tags

| Category | Tag | Description |
| :--- | :--- | :--- |
| **Tips** | `linear_with_tip` | A path ending in a dead-end node. Usually represents sequencing errors or low-coverage regions. |
| **Circular** | `circular_isolated` | A standalone perfect loop with no entry/exit points in the graph. Often represents plasmids or organelle genomes. |
| | `circular_fork_hub_residual` | A circular path extracted during the residual processing of branching nodes. |
| | `circular_with_bubble` | A circular structure that incorporates a resolved variant (bubble). |
| | `circular_self_loop_hub` | A path that starts and ends at the same hub node without any intermediate branches. |
| **Lariat** | `lariat_simple` | A 'lasso' topology: a linear 'tail' leading into a cycle. |
| | `lariat_branch_at_hub` | A lariat structure where the junction point is a recognized graph hub. |
| | `lariat_fork_hub_residual` | A lariat path identified during residual hub cleanup. |
| **Linear** | `linear_simple` | A standard linear path with terminal ends and no internal branching. |
| | `linear_bubble_resolved_ending` | A branch of a resolved bubble that terminates exactly at the convergence node. |
| | `linear_bubble_resolved_extended` | A branch of a resolved bubble that continues past the convergence point. |
| | `linear_through_hub` | A linear path that traverses through one or more branching points. |



---

## 2. Dynamic Suffixes (Connectivity Map)

Most tags (except for isolated circles and tips) are appended with suffixes that describe their start and end conditions relative to the graph's complexity.

| Suffix | Logic | Meaning |
| :--- | :--- | :--- |
| `_bubble_starting` | `start_node == bubble_source` | The contig originates exactly at the bifurcation point of a resolved variant. |
| `_hub_starting` | `start_node.is_hub()` | The contig originates from a generic branching node (In-degree > 1 or Out-degree > 1). |
| `_hub_ending` | `final_node.is_hub()` | The contig terminates at a branching node (convergence or new divergence). |

---

## 3. Unknown/Ambiguous Tags
These tags are used for paths that don't fit the strict definitions above, often appearing during the final "catch-all" search for remaining graph components.
* `unknown_circular_like`
* `unknown_lariat_like`
* `unknown_linear_like_from_catch_all`

---

## Example Header Analysis
`>contig_42|len=2500|tag=linear_bubble_resolved_ending_bubble_starting_hub_ending|cov=30.0`

**Interpretation:**
- **Structure:** It is a linear branch of a resolved bubble.
- **Start:** It began exactly at the bubble's source node (`_bubble_starting`).
- **End:** It finished at the bubble's destination node, which is also a hub (`_hub_ending`).
- **Significance:** High-confidence allelic variant or paralogous segment.