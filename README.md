# README WORK IN PROGRESS
# De Bruijn Graph Assembler & Contig Classifier

This tool implements a genome assembler based on the De Bruijn Graph formalism. It accepts k-mer inputs to construct the graph structure and performs a topology-aware assembly process.

Unlike standard greedy assemblers, this implementation explicitly identifies and handles graph complexities such as sequencing artifacts (tips), allelic variants (bubbles), and isolated components. The assembly engine traces paths through the graph and assigns a specific topological class to each generated contig, providing granular insights into the graph structure directly in the output metadata.
## Key Features

- **Structural Variant Resolution:** Automated handling of *bubbles* (allelic variants/paralogs) and *tips* (sequencing artifacts) during path extension.
- **Advanced Topological Tagging:** Every contig is classified based on its connectivity and shape within the graph (e.g., Linear, Circular, Lariat).
- **Metadata-Rich FASTA Output:** Generates FASTA headers containing median k-mer coverage, length, topological tags and necessary structural information to better understand the origin of the topological complications.




# Contig  Tagging

Assembled contigs are classified according to their origin and topology. Every sequence is tagged based on its **shape** in the graph and the particular structures (like bubbles or tips) it contains.

## 1. Primary Topologies
First, the contig is put into one of three main families based on its start and end points.

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
If it is not a tip, the code looks at the shape and where it came from:

* **Circular**:
  * `circular_isolated`: A perfect loop found at the very end (no junctions).
  * `circular_from_hub`: A loop that started from a junction (hub).
* **Lariat**:
  * `lariat_hub_branching`: The lasso loop closes at a junction (hub).
  * `lariat`: The lasso loop closes at a simple node.
* **Linear**:
  * `linear`: The default tag for open-ended paths.

### C. Added Suffixes
Unless it is a tip, we add extra info to the base tag to describe the journey:
* **`_with_bubble`**: The contig solved a genetic variant (bubble).
* **`_through_hub`**: The path passed through at least one junction (hub).
* **`_simple`**: A clean path with no bubbles and no junctions.

---

## 3. FASTA Header Format

Each sequence is saved with a header containing metadata separated by pipes (`|`):

`>contig_ID|len=BP|tag=TAG|cov=MEDIAN|assembly_phase=PHASE|...`

### Extra Metadata (Conditional)
Additional info appears only if it fits the structure:
* **`start/end`**: Tells you if the end is a `hub`, a `bubble_start/end`, or a `tip_start`.
* **`branch`**: For lariats, shows the k-mer where the tail meets the loop.
* **`bubble_ends`**: Shows the start and end k-mers of a solved variant.
* **`tip_from`**: Shows the junction where a tip originally branched off.
* **`hub_origin`**: The main junction where the search started.

---

## 4. Quick Technical Notes
1. **Median Coverage**: The `cov` value uses the **median** weight. This is better than the average because it isn't ruined by small repetitive peaks.
2. **The Catch-all Rule**: `circular_isolated` contigs are almost always `_simple`. If they had a bubble, they would have a junction and would have been found earlier in the process.
3. **Adding Contigs**: You can use the `+` operator to merge two contigs. The new contig will combine all flags (like `bubble` or `tip`) and update the list of visited junctions.