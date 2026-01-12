# README WORK IN PROGRESS
# De Bruijn Graph Assembler & Contig Classifier

This tool implements a genome assembler based on the De Bruijn Graph formalism. It accepts k-mer inputs to construct the graph structure and performs a topology-aware assembly process.

This implementation can explicitly identify and handle graph complexities such as sequencing artifacts (tips or bubbles), allelic variants (bubbles), and isolated components. The assembly engine traces paths through the graph and assigns a specific topological class to each generated contig. This provides granular insights into the graph structure directly in the output metadata, enabling downstream analysis using comprehensive parsing-ready FASTA headers.
## Key Features

- **Explicit Stats calculation** At any moment you can ask for the current graph stats, that give you insight on the actual number of *edges*, *nodes* (of which *hubs* and *receiving hubs*), *average coverage* and number of *bubbles and tips* kept after processing
- **Eulerian Path Detection** Explicitly detects the presence of a single *eulerian path* (circular or linear), fundamental to understand if the cleaning graph procedure managed to fully simplify the graph (if possible)
- **Structural Variant Resolution:** Automated handling of *bubbles* (allelic variants/sequencing artifacts) and *tips* (sequencing artifacts) during path extension. Bubbbles and tips can be removed or kept if considered biological thruths, through user-defined thresholds
- **fork_in_hub option** When enabled, this parameter allows to explore one level of depthness of the graph, allowing to resolve bubble structures or to walk through hubs keeping the 5' part of the contig as common between different complex contigs. This, however, poses the risk of creating chimeric assemblies and can overstimate the n50 calculation
- **Advanced Topological Tagging:** Every contig is classified based on its connectivity and shape within the graph (e.g., *Linear, Circular, Lariat*), the assembly phase where it generated and where it started/ends
- **Metadata-Rich FASTA Output:** Generates FASTA headers containing median k-mer coverage, length, topological tags and necessary structural information to better understand the origin of the topological complications that interrupted the assembly.

## USAGE

This assembler can be used in two ways: as a python package and using a pre-defined highly customizable pipeline function.
To start, you can simply clone this repository.




# Contig  Tagging

Assembled contigs are classified according to their origin and topology. Every sequence is tagged based on its **shape** in the graph and the particular structures (like bubbles or tips) it contains.

## 1. Primary Topologies
First, the contig is put into one of three main families (`main_tag=` in the FASTA header) based on its start and end points.

| Family | Logic | Description |
| :--- | :--- | :--- |
| **`linear`** | Start and end are different. | A standard piece of DNA. |
| **`circular`** | Start and end are the same node. | A complete circular element (like a plasmid). |
| **`lariat`** | The end node was already visited. | A "lasso" shape (a tail leading into a loop). |



---

## 2. How Tags are Assigned (`get_tag`)

The `get_tag()` method builds the final name using a strict order of importance:

### A. Top Priority: Tips
If the contig is a "tip" (a dead-end), it is immediately tagged as:
* **`linear_with_tip`**
> *Note: If it’s a tip, the analysis stops here. No other suffixes are added.*

### B. Base Tag (If not a Tip)
If it is not a tip, the code looks at the shape and gives additional information for specific circular or lariat shapes:

* **Circular**:
  * `circular_isolated`: A perfect loop found at the very end (no junctions). It starts from a node that is no hub and has one edge incoming. It must have been generated during the `catch_all` phase
  * `circular_from_hub`: A loop that started from a junction (hub). 
* **Lariat**:
  * `lariat_hub_branching`: The lasso loop closes at a junction (hub).
  * `lariat`: The lasso loop closes at a simple node.
* **Linear**:
  * `linear`: The default tag for open-ended paths.

### C. Added Suffixes
Unless it is a tip, we add extra info to the base tag to describe the walk that created the contig:
* **`_with_bubble`**: The contig solved a bubble.
* **`_through_hub`**: The path passed through at least one junction (hub).
* **`_simple`**: A clean path with no bubbles and no junctions.

> * Note: `_with_bubble` and `_through_hub` are compatible for bubble resolutions that result in a hub as the final bubble node.

---

## 3. FASTA Header Format

Each sequence is saved with a header containing metadata separated by pipes (`|`):

`>contig_ID|len=BP|main_tag=MAIN_TAG|tag=TAG|cov=MEDIAN|assembly_phase=PHASE|EXTRA_METADATA`

### Median Coverage:

The `cov` value uses the **median** weight. This is better than the average because it isn't ruined by small repetitive peaks.

### Assembly Phase
The phase of the assembly process in which the contig generated:
* **`standard`**: The initial contig assembly phase, starting only from nodes with no incoming edges if `fork_in_hub` is set to True, otherwise also all hubs are considered as initial starts.
* **`residual`**: Possible only when `fork_in_hub` is True. It is assigned to all contigs starting from hubs where the initial standard search stopped because limited to one level of depthness
* **`from_catch_all`**: assigned to all contigs generated in the final phase of assembly to detect contigs not starting neither from nodes with no incoming edges nor hubs. Should be assigned only to `circular_isolated` tagged


### Extra Metadata (Conditional)
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

## 4. Quick Technical Notes

**Adding Contigs**: You can use the `+` operator to merge two contigs. The new contig will combine all flags (like `bubble` or `tip`) and update the list of visited junctions. However, you can not add contigs originated from different assembly phases.
