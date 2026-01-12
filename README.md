# README WORK IN PROGRESS
# De Bruijn Graph Assembler & Contig Classifier

This tool implements a genome assembler based on the De Bruijn Graph formalism. It accepts k-mer inputs to construct the graph structure and performs a topology-aware assembly process.

This implementation can explicitly identify and handle graph complexities such as sequencing artifacts (tips or bubbles), allelic variants (bubbles), and isolated components. The assembly engine traces paths through the graph and assigns a specific topological class to each generated contig. This provides granular insights into the graph structure directly in the output metadata, enabling downstream analysis using comprehensive parsing-ready FASTA headers.

---
## Key Features

- **Explicit Stats calculation** At any moment you can ask for the current graph stats, that give you insight on the actual number of *edges*, *nodes* (including *hubs* and *receiving hubs*), *average coverage* and number of *bubbles and tips* kept after processing
- **Eulerian Path Detection** Explicitly detects the presence of a single *eulerian path* (circular or linear), fundamental to understand if the cleaning graph procedure managed to fully simplify the graph (if possible)
- **Structural Variant Resolution:** Automated handling of *bubbles* (allelic variants/sequencing artifacts) and *tips* (sequencing artifacts) during path extension. Bubbles and tips can be removed or kept if considered biological variants, through user-defined thresholds
- **fork_in_hub option** When enabled, this parameter allows to explore one level of depth of the graph, allowing to resolve bubble structures or to walk through hubs keeping the 5' part of the contig as common between different complex contigs. This, however, poses the risk of creating chimeric assemblies and can overestimate the n50 calculation
- **Advanced Topological Tagging:** Every contig is classified based on its connectivity and shape within the graph (e.g., *Linear, Circular, Lariat*), the assembly phase where it was generated and where it started/ended
- **Metadata-Rich FASTA Output:** Generates FASTA headers containing median k-mer coverage, length, topological tags and necessary structural information to better understand the origin of the topological complications that interrupted the assembly.

---

## 📚 Key Terms

### Graph Components
- Node: A (k-1)-mer sequence. Nodes are connected by edges.
- Edge: A k-mer with a weight (coverage) indicating how many times it appears in the reads. It connects two nodes starting from the node containing its -1 prefix and arriving at the node containing its -1 suffix
- Hub: A node with multiple outgoing paths (branching point).
- Receiving Hub: A node where multiple paths converge.

### Graph Particular Structures
- Bubble: Divergence + reconvergence pattern. Represents variants or sequencing errors.
- Tip: Short dead-end branch, usually from sequencing errors at read ends. We distinguish type I (outgoing) tips as tips that start from a hub and result in a dead end and type II (incoming) for those tips that start from a node with no incoming edges and shortly after converge in a main structure.
- Island: Isolated graph component, disconnected from the main structure.


> **Note**: All those structures are detected only when their length is **shorter than 2k-1 bases**. If longer, it would be hard to consider them as errors and more probable to consider them as biologic variants. When detected, they can be removed or kept in the graph via user defined thresholds.



### Contig Shapes
- Linear: Standard path with distinct start/end.
- Circular: Closed loop (start = end). Often plasmids.
- Lariat: "Lasso" shape - linear tail into a cycle.

### Assembly Terms
- Coverage: Number of reads supporting a k-mer. Higher = more confident.
- Eulerian Path: Path visiting every edge exactly once. Perfect assemblability.
- N50: Assembly quality metric. Length where 50% of bases are in contigs ≥ that length.

---

## USAGE

The assembler is designed for both automated execution and manual graph manipulation. You can run the entire pipeline at once or control each step through the Dbg class.

### 1. Automated Pipeline (Quick Start)
The run_assembly_pipeline function is the fastest way to assemble your data. It handles everything from FASTQ parsing to N50 calculation. Example usage is below.

```python
from pipeline_function import run_assembly_pipeline

dbg, lengths, n50 = run_assembly_pipeline(
    k=31,
    file_name="reads_R1.fq",
    file_type="fastq",
    paired=True,
    file_nameR2="reads_R2.fq",
    quality_check=True,
    prune_threshold=3,          # Removes low-frequency k-mers (noise)
    tip_cov_threshold=4,        # Coverage threshold for tip removal
    bubble_cov_threshold=4,     # Coverage threshold for bubble resolution
    fork_in_hub=True,           # Enables deeper path exploration
    output_fasta="assembly_output.fasta"
) 
```

You can find full documentation for this function in its file.


### 2. Manual Control (Step-by-Step)

For complete control over the assembly process, you can execute each step individually using the Dbg class. See `main_pipeline.py` for a complete working example.

Basic Workflow:

1. **Initialize** → Create the Dbg object with your k-mer size and input file(s)
2. **Build Graph** → Parse reads and construct the De Bruijn Graph
3. **Check Stats** → Inspect the initial graph structure (nodes, edges, hubs)
4. **Clean Graph** → Remove noise and resolve structural issues:
   - **prune_graph()** → Remove low-coverage k-mers (sequencing errors)
   - **process_tips()** → Remove dead-end branches (read-end errors)
   - **process_bubbles()** → Resolve divergence-reconvergence patterns (variants/errors)
   - **process_short_islands()** → Remove small isolated components 
5. **Check Stats Again** → Verify the cleaning improved the graph
6. **Check if graph is eulerian** → Verify if the cleaning put the graph in eulerian conditions: perfect assembly expected
7. **Assemble** → Generate contigs by traversing the cleaned graph
8. **Write Output** → Save contigs to FASTA with metadata
9. **Calculate N50** → Assess assembly quality

**Key Differences from Automated Pipeline**:
- You can inspect graph statistics between steps to guide parameter choices
- You can skip or repeat cleaning steps (e.g., run process_bubbles() twice with different thresholds)
- You can access intermediate results for debugging or analysis (e.g. create your custom stats)
- You can customize the order of operations
- You can perform different assemblies before and after some cleaning steps, or both using and not using `fork_in_hub` option.

**When to Use Manual Control**:
- Optimizing parameters
- Debugging assembly issues
- Running custom analyses on the graph structure

---


## ⚙️ Parameter Guide

### Graph Building

**`k` (k-mer size)**
- **Lower k (15-25)5)**: More connections, better for low coverage or diverse data. Risk: more branches
- Higher k (31-55)5)**: Fewer spurious connections, better for high coverage. Risk: fragmentation
> Advice: Use odd k, start with 31.

**`quality_threshold` (Phred score)**
- **Higher (30+)**: Only high-quality reads. More accurate but less 
- **Lower (20-25)**: More data but more noise. Use with higher pruning.
> **NOTE**: This assembler **is not** a trimmer. The quality check logic is very brutal by discarding whole reads if their median is not above the threshold. Best advice is using this assembler on already trimmed data with more reasonable logic, or using 30 or less as threshold.

**`reverse_r2`**
- **True**: Considers the reads of the R2 file as already reverse complemented by the assembler. Makes sense if you have data you want to use as stranded and you know this information.
- **False** (default): considers the R2 file reads as on the opposite strand of the R1 file. Is perfectly ok if you want to use both strands however. See below.

**`use_both_strands`**
- **True** (default): Considers both DNA strands. Recommended unless you know reads are strand specific and you need this feature.
- **False**: Only forward strand. Rare use case.


### Graph Cleaning

**`prune_graph(threshold)`**
- Removes edges with coverage < threshold.
- **threshold=2**: Removes singleton k-mers (likely errors)
- **threshold=3-5**: If you know your data can be a bit noisy
- **Higher values**: Risk losing real low-coverage regions.

**`process_tips(cov_threshold)`**
- Removes short dead-end branches (tips) with low coverage.
- **cov_threshold=5**: Standard. Keeps tips with ≥5X coverage, that have probably some biological basis.
- **Higher values**: More aggressive, can potentially remove biological information but simplify the graph.
- **Lower values**: Keep more potential biology, more noise.

**`process_bubbles(cov_threshold, ratio_threshold)`**
- Resolves divergence-convergence patterns (bubbles).
- **cov_threshold**: Minimum coverage to keep a bubble path.
  - 5: Keep paths with ≥5X (likely real variants). Default value.
  - Higher: Only high-confidence variants.
  - Lower: More permissive, can keep some noise.
- **ratio_threshold**: Minimum ratio to the best path.
  - 0.1: Keep paths with ≥10% of best path coverage.
  - 0.2: Stricter, fewer paths kept.
  - 0.05: More permissive.
> **General rule**: Lower values keep more variants (good for biological variability); higher values simplify the graph but can lose information.

**`process_short_islands()`**
- Removes isolated components shorter than 2k+1.
- No parameters - automatically removes likely contamination

**`fork_in_hub`**
- **False** (default): Conservative. Stops at hubs, each path becomes separate contig.
  - **Use when**: You want accurate assemblies that cannot risk having chimeric contigs.
  - **Best for**: Circular genome detection.
  
- **True**: Explores one level beyond hubs. Resolves all bubble paths and different hub paths resulting in contigs that share a common head 5' portion.
  - **Use when**: Complex region with variants you want to investigate.
  - **Warning**: Can create chimeric assemblies. N50 can be overestimated.
  - **Best for**: Maximizing contig length, resolving structural variants

**`write_contigs(min_len)`**
 
 Only saves contigs having length >= min_len.
- **min_len=200**: Remove very short fragments.
- **min_len=500-1000**: Focuses on very relevant contigs.
- **min_len=0**: Save everything (useful for debugging).

---

## Contig Tagging

Assembled contigs are classified according to their origin and topology. Every sequence is tagged based on its **shape** in the graph and the particular structures (like bubbles or tips) it contains.

### 1. Primary Topologies (main_tag)
First, the contig is put into one of three main families (`main_tag=` in the FASTA header) based on its start and end points.

| Family | Logic | Description |
| :--- | :--- | :--- |
| **`linear`** | Start and end are different. | A standard piece of DNA. |
| **`circular`** | Start and end are the same node. | A complete circular element (like a plasmid). |
| **`lariat`** | The end node was already visited. | A "lasso" shape (a tail leading into a loop). |



---

### 2. How Tags are Assigned (`get_tag`)

The `get_tag()` method builds the final name using a strict order of importance:

#### A. Top Priority: Tips
If the contig is a "tip" (a dead-end), it is immediately tagged as:
* **`linear_with_tip`**
> **Note**: If it’s a tip, the analysis stops here. No other suffixes are added.*

#### B. Base Tag (If not a Tip)
If it is not a tip, the code looks at the shape and gives additional information for specific circular or lariat shapes:

* **Circular**:
  * `circular_isolated`: A perfect loop found at the very end (no junctions). It starts from a node that is no hub and has one edge incoming. It must have been generated during the `catch_all` phase
  * `circular_from_hub`: A loop that started from a junction (hub). 
* **Lariat**:
  * `lariat_hub_branching`: The lasso loop closes at a junction (hub).
  * `lariat`: The lasso loop closes at a simple node.
* **Linear**:
  * `linear`: The default tag for open-ended paths.

#### C. Added Suffixes
Unless it is a tip, we add extra info to the base tag to describe the walk that created the contig:
* **`_with_bubble`**: The contig solved a bubble.
* **`_through_hub`**: The path passed through at least one junction (hub).
* **`_simple`**: A clean path with no bubbles and no junctions.

> **Note**: `_with_bubble` and `_through_hub` are compatible for bubble resolutions that result in a hub as the final bubble node.

---

### 3. FASTA Header Format

Each sequence is saved with a header containing metadata separated by pipes (`|`):

`>contig_ID|len=BP|main_tag=MAIN_TAG|tag=TAG|cov=MEDIAN|assembly_phase=PHASE|EXTRA_METADATA`

#### Median Coverage:

The `cov` value uses the **median** weight. This is better than the average because it isn't ruined by small repetitive peaks.

#### Assembly Phase
The phase of the assembly process in which the contig was generated:
* **`standard`**: The initial contig assembly phase, starting only from nodes with no incoming edges, if `fork_in_hub` is set to True, otherwise also all hubs are considered as initial starts.
* **`residual`**: Possible only when `fork_in_hub` is True. It is assigned to all contigs starting from hubs where the initial standard search stopped because limited to one level of depth
* **`from_catch_all`**: assigned to all contigs generated in the final phase of assembly to detect contigs starting neither from nodes with no incoming edges nor from hubs. Should be assigned only to contigs tagged as `circular_isolated`


#### Extra Metadata (Conditional)
Additional info appears only if it fits the structure:
* **`start/end`**: Tells you if the end is a `hub`, a `bubble_start/end`, or a `tip_start`.
* **`branch`**: For lariats, shows the k-mer where the tail meets the loop.
* **`bubble_ends`**: Shows the start and end k-mers of a solved variant.
* **`tip_from`**: Shows the junction where a tip originally branched off.
* **`hub_seq`**: The main junction where the search forked (when `fork_in_hub = True`).


**Important note on circular contigs**:

Some tagging combinations should not be possible and represent an immediate warning when found: **circular contigs can never be `_through_hub`** because they either start from hub (in the residual assembly phase or when `fork_in_hub` is False, and in both conditions the assembly stops at hubs) or are caught in the catch all final phase when no hubs are present. 

They can however be `_with_bubble` when the hub they start from is directly a bubble start node.

Finally, `circular_isolated` contigs should always be`_simple`. If they had a bubble, they would have a junction and would have been found earlier in the process.

Therefore, if you know that you are specifically searching for circular contigs you are advised to use `fork_in_hub=False` allowing the detection of circular starting from primary and not nested hubs.


---

### 4. Quick Technical Notes

**Adding Contigs**: You can use the `+` operator to merge two contigs. The new contig will combine all flags (like `bubble` or `tip`) and update the list of visited junctions. However, you cannot add contigs originated from different assembly phases.
