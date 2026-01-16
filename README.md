# De Bruijn Graph Assembler & Contig Classifier

This tool implements a genome assembler based on the De Bruijn Graph strategy. It accepts k-mer inputs (calculated from FASTQ or FASTA parsed files) to construct the graph structure and performs a topology-aware assembly process.

This implementation can explicitly identify and handle graph complexities such as sequencing artifacts (tips or bubbles), biological variants (bubbles), and isolated components. The assembly engine traces paths through the graph and assigns a specific topological tag to each generated contig.

This provides granular insights into the graph structure directly in the output metadata, enabling downstream analysis using comprehensive parsing-ready FASTA headers.

---
## 📌 Table of Contents
- [Key Features](#key-features)
- [📚 Key Terms](#-key-terms)
  - [Graph Components](#graph-components)
  - [Graph Particular Structures](#graph-particular-structures)
  - [Contig Shapes](#contig-shapes)
  - [Assembly Terms](#assembly-terms)
- [Installation](#installation)
- [📖 Technical Documentation](#-technical-documentation)
- [USAGE](#usage)
  - [1. Automated Pipeline (Quick Start)](#1-automated-pipeline-quick-start)
  - [2. Manual Control (Step-by-Step)](#2-manual-control-step-by-step)
- [⚙️ Parameter Guide](#%EF%B8%8F-parameter-guide)
- [Contig Tagging](#contig-tagging)
  - [1. Primary Topologies](#1-primary-topologies-main_tag)
  - [2. How Tags are Assigned](#2-how-tags-are-assigned-get_tag)
  - [3. FASTA Header Format](#3-fasta-header-format)
- [Areas for Improvement](#areas-for-improvement)
- [👥 Credits](#-credits)



---
## Key Features

- **Explicit Stats calculation** At any moment you can ask for the current graph stats, that give you insight on the actual number of *edges*, *nodes* (including *hubs* and *receiving hubs*), *average coverage* and number of *bubbles and tips* kept after processing
- **Eulerian Path Detection** Explicitly detects the presence of a single *eulerian path* (circular or linear), fundamental to understand if the cleaning graph procedure managed to fully simplify the graph (if possible)
- **Structural Variant Resolution:** Automated handling of *bubbles* (biological variants/sequencing artifacts) and *tips* (sequencing artifacts) during path extension. Bubbles and tips can be removed or kept if considered biological variants, through user-defined thresholds
- **fork_in_hub option** When enabled, this parameter allows to explore one level of depth in the graph, allowing to resolve bubble structures or to walk through hubs keeping the 5' part of the contig as common between different complex contigs. This, however, poses the risk of creating chimeric assemblies and can overestimate the n50 calculation
- **Advanced Topological Tagging:** Every contig is classified based on its connectivity and shape within the graph (e.g., *Linear, Circular, Lariat*), the assembly phase where it was generated and where it started/ended
- **Metadata-Rich FASTA Output:** Generates FASTA headers containing median k-mer coverage, length, topological tags and necessary structural information to better understand the origin of the topological complications that interrupted the assembly.

---

## 📚 Key Terms

### Graph Components
- **Node**: A (k-1)-mer sequence. Nodes are connected by edges.
- **Edge**: A k-mer with a weight (coverage) indicating how many times it appears in the reads. It connects two nodes: the source node (containing its first k-1 bases, or prefix) and the destination node (containing its last k-1 base, or suffix)
- **Hub**: A node with multiple outgoing paths (branching point).
- **Receiving Hub**: A node where multiple paths converge.

### Graph Particular Structures
- **Bubble**: Divergence + reconvergence pattern. Represents variants or sequencing errors.
- **Tip**: Short dead-end branch, usually from sequencing errors at read ends. We distinguish type I (outgoing) tips as tips that start from a hub and result in a dead end and type II (incoming) for those tips that start from a node with no incoming edges and shortly after converge in a main structure.
- **Island**: Isolated graph component, disconnected from the main structure.


> **Note**: All those structures are detected only when their length is **shorter than 2k+1 bases**. If longer, it would be hard to consider them as errors and more probable to consider them as biological variants. When detected, they can be removed or kept in the graph via user defined thresholds.



### Contig Shapes
- **Linear**: Standard path with distinct start/end.
- **Circular**: Closed loop (start = end). 
- **Lariat**: "Lasso" shape - linear tail into a cycle.

### Assembly Terms
- **Coverage**: Number of reads supporting a k-mer. Higher = more confident.
- **Eulerian Path**: Path visiting every edge exactly once. Perfect assemblability.
- **N50**: Assembly quality metric. Length where 50% of bases are in contigs ≥ that length.

---

## Installation

### Requirements
- Python 3.9 or higher
- No external dependencies (uses only Python standard library)

### Setup
```bash
git clone git@github.com:alfalconetti/DBG_Assembly_Bioinfo25_Civilwar.git

cd DBG_Assembly_Bioinfo25_Civilwar

```

All modules must be imported from `src` folder.

--- 

## 📖 Technical Documentation
All documentation is available via Sphinx using a browser based interface. 

To generate and view the documentation:
1. Ensure you have the theme installed: `pip install sphinx-rtd-theme`
2. Run the build command:
   ```bash
   cd docs
   make html
	```

​Open `docs/_build/html/index.html` in your browser.

You can also use the `documentation.pdf` in the root repo folder but it may not be up to date

​This will provide a searchable and interconnected interface for all classes and their methods.

---

## USAGE

The assembler is designed for both automated execution and manual graph manipulation. You can run the entire pipeline at once or control each step through the Dbg class.
The starting point can be FASTA or FASTQ files.
Only FASTQ files can be paired.
The assembler also is able to read compressed .gz files.

### 1. Automated Pipeline (Quick Start)
The run_assembly_pipeline function is the fastest way to assemble your data. It handles everything from FASTQ parsing to N50 calculation. Example usage is below.

```python
from src.pipeline_function import run_assembly_pipeline

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

You can find full documentation for this function as described above

A working example is in the `main_pipeline_function.py` file in the root repo directory

A Command line interface for this function could be developed to ease its use in bash workflows.


### 2. Manual Control (Step-by-Step)

For complete control over the assembly process, you can execute each step individually using the Dbg class. See `main_pipeline.py` for a complete step by step working example.

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
   - **`is_eulerian_linear()`** 
   - **`is_eulerian_circular()`**
7. **Assemble** → Generate contigs by traversing the cleaned graph
8. **Write Output** → Save contigs to FASTA with metadata
9. **Calculate N50** → Assess assembly quality

**Key Differences from Automated Pipeline**:
- You can inspect graph statistics between steps to guide parameter choices
- You can skip or repeat cleaning steps (e.g., run `process_bubbles()` twice with different thresholds)
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
- **Lower k (15-25)**: More connections, better for low coverage or diverse data. Risk: more branches
- **Higher k (31-55)**: Fewer spurious connections, better for high coverage. Risk: fragmentation
> Advice: Use odd k, start with 31.

**`quality_threshold` (Phred score)**
- **Higher (30+)**: Only high-quality reads. More accurate but less information
- **Lower (20-25)**: More data but more noise. Use with higher pruning.
> **NOTE**: This assembler **is not** a trimmer. The quality check logic is very brutal by discarding whole reads if their median is not above the threshold. Best advice is using this assembler on already trimmed data with more reasonable trimming logic, or using 30 or less as threshold.

**`reverse_r2`**
- **True**: Considers the reads of the R2 file as already reverse complemented during the sequencing. Makes sense if you have data you want to use as stranded and you know this information.
- **False** (default): considers the R2 file reads as on the opposite strand of the R1 file. Default is perfectly ok if you want to use both strands however. See below.

**`use_both_strands`**
- **True** (default): Considers both DNA strands. Recommended unless you know reads are strand specific and you need this feature.
- **False**: Only forward strand. Rare use case.

> **Note**: If `use_both_strands=True` it overrides `reverse_r2` because both strands are used anyway.


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
- No parameters - automatically removes likely errors

**`fork_in_hub`**
- **False** (default): Conservative. Stops at hubs, each path becomes separate contig.
  - **Use when**: You want accurate assemblies that cannot risk having chimeric contigs.
  - **Best for**: Circular genome detection.
  - **Result**: More contigs covering independent assembly regions
  
- **True**: Explores one level of depth beyond hubs. When it encounters an hub, it saves the contig up to that point. Then, it resolves all bubble paths and different hub paths and sums them to the saved contig, resulting in contigs that share a common head 5' prefix. 
  - **Use when**: Complex region with variants you want to investigate.
  - **Warning**: Can create chimeric assemblies. N50 can be overestimated.
  - **Best for**: Maximizing contig length, resolving structural variants


> Technical detail: "One level" means the assembler explores only the immediate children paths of the hub, not their children. This prevents exponential path explosion while still resolving local complexity. 
The most complicated scenario it handles (very rare) is the presence of an hub in the final node of a bubble that creates a number of contigs equal to the product of the bubble paths and the subsequent hub paths.

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


### 3. FASTA Header Format

Each sequence is saved with a header containing metadata separated by pipes (`|`):

`>contig_ID|len=BP|main_tag=MAIN_TAG|tag=TAG|cov=MEDIAN|assembly_phase=PHASE|EXTRA_METADATA`

#### Median Coverage:

The `cov` value uses the **median** weight. This is better than the average because it isn't ruined by small repetitive peaks.

#### Assembly Phase
The phase of the assembly process in which the contig was generated:
* **`standard`**: The initial contig assembly phase, starting only from nodes with no incoming edges, if `fork_in_hub` is set to True, otherwise also all hubs are considered as initial starts.
* **`residual`**: Possible only when `fork_in_hub` is True. It is assigned to all contigs starting from hubs where the initial standard search stopped because limited to one level of depth (secondary hubs)
* **`from_catch_all`**: assigned to all contigs generated in the final phase of assembly to detect contigs starting neither from nodes with no incoming edges nor from hubs. Should be assigned only to contigs tagged as `circular_isolated`


#### Extra Metadata (Conditional)
Additional info appears only if it fits the structure:
* **`start/end`**: Tells you if the end is a `hub`, a `bubble_start/end`, or a `tip_start`.
* **`branch`**: For lariats, shows the k-mer where the tail meets the loop.
* **`bubble_ends`**: Shows the start and end k-mers of a solved variant.
* **`tip_from`**: Shows the junction where a tip originally branched off.
* **`hub_seq`**: The main junction where the search forked (when `fork_in_hub = True`).


**Important note on circular contigs**:

Some tagging combinations should not be possible and represent an immediate warning when found: **circular contigs can never be `_through_hub`** because they either start from hub (in the residual assembly phase or when `fork_in_hub` is False, and in both conditions the assembly stops as soon as it finds hubs) or are caught in the catch all final phase when no hubs are present. 

They can however be `_with_bubble` when the hub they start from is directly a bubble start node.

Finally, `circular_isolated` contigs should always be`_simple`. If they had a bubble, they would have a junction and would have been found earlier in the process.

> **ADVICE**: if you know that you are specifically searching for circular contigs you are advised to use `fork_in_hub=False` allowing the detection of circular starting from primary and not nested hubs.


### 4. Quick Technical Notes

**Adding Contigs**: You can use the `+` operator to merge two contigs. The new contig will combine all flags (like `bubble` or `tip`) and update the list of visited junctions. However, you cannot add contigs originated from different assembly phases. 

---

##  Areas for Improvement

While this assembler provides a robust framework for De Bruijn Graph construction and topological analysis, several enhancements could improve its performance and scalability:

### 1. Parallelization in Graph Construction
Currently, the graph is built sequentially. Implementing a multi-threaded approach for k-mer counting and edge addition would significantly reduce the execution time, especially for large multi-gigabyte FASTQ files. This has not been done yet because the Python GIL locks from accessing the same memory area (the nodes and edges dictionaries), this would force to create separate dictionaries to merge and the merging process may be quite computationally intensive.

### 2. Scaffolding with Paired-End Reads
The current engine focuses on contig assembly through graph linear walk that stops in hubs or at nodes with no outgoing edges. Integrating paired-end information during or after the assembly process would allow the join different contigs, resolving repeats and creating longer, more continuous **scaffolds**. This is the case for more complex genomes that may result in an highly fragmented assembly (low N50)

### 3. Graph Compaction with Supernodes
To optimize memory usage and simplify traversal, the graph could be compacted into a **Compressed De Bruijn Graph**. By merging non-branching paths (chains of nodes with in-degree and out-degree of 1, that is almost equivalent to our non-forked contigs) into single **supernodes**, the number of elements to manage would drastically decrease without losing genomic information. This is totally useless for simple small genomes but is of uttermost importance if this tool were to be used in big complex genomes, in order to reduce the impact of the tool on memory.



### 4. Command Line Interface (CLI)
Developing a dedicated CLI (e.g. using `argparse`) would make the tool more accessible and easier to integrate into automated bioinformatics pipelines, allowing users to specify parameters, input files, and output paths directly from the terminal, offering another highly customizable usage option.


---

## 👥 Credits

This project was developed as part of the Programmazione per la Bioinformatica course (2025/26) in University of Bari.

* *Developers*
  * De Lorenzis Stefano - [stedelorenzis]
  * Falconetti Andrea Luca - [alfalconetti]

