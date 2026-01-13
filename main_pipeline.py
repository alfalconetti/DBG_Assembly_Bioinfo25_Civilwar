import time
import sys
from src.dbg import Dbg


def main():
    """
        The main execution entry point for the De Bruijn Graph Assembler Pipeline.

        This function coordinates the entire lifecycle of the assembly process,
        from reading raw sequencing data to producing the final FASTA output.

        The pipeline is structured in 5 sequential stages:
        1.  Graph Construction: Parses FASTQ files and builds the initial De Bruijn Graph.
        2.  Pruning: Removes k-mers with low coverage (likely sequencing errors) to clean the graph.
        3.  Topological Cleaning: Resolves complex graph structures:
                - Tips: Dead ends caused by errors at the end of reads.
                - Bubbles: Diverging paths caused by SNPs or errors in the middle of reads.
                - Islands: Disconnected small components.
        4.  Assembly: Traverses the cleaned graph to generate contigs (sequences).
        5.  Output & Stats: Writes the results to a file and calculates quality metrics (N50).

        Configuration:
        File paths and parameters (k-mer size, thresholds) are currently defined
        as constants at the beginning of this function.

        :raises SystemExit: If a critical error occurs during graph construction (e.g., file not found).
    """
    FILE_R1 = "test_data/Illumina_reads1.fq"
    FILE_R2 = "test_data/Illumina_reads2.fq"
    K_MER_SIZE = 31
    OUTPUT_FILE = "assembly_output.fasta"

    print(f"Starting De Bruijn Graph Assembler Pipeline...")
    start_global = time.time()


    # 1. Graph construction
    print(f"\n[1] GRAPH CONSTRUCTION (k={K_MER_SIZE})")
    print(f"Dataset: Paired-End ({FILE_R1}, {FILE_R2})")

    t0 = time.time()
    try:
        grafo = Dbg(K_MER_SIZE, FILE_R1, "fastq", paired=True, file_nameR2=FILE_R2)
        print("> Parsing reads and building nodes...")
        grafo.build_graph(quality_check=True, quality_threshold=20)
    except Exception as e:
        print(f"Critical Error during construction: {e}")
        sys.exit(1)

    print(f"> Time elapsed: {time.time() - t0:.2f}s")
    grafo.get_graph_stats()

    # 2. Pruning
    print(f"\n[2] PRUNING (Noise Reduction)")
    print("> Removing k-mers with low frequency (threshold=2)...")

    t0 = time.time()
    grafo.prune_graph(threshold=3)
    print(f"> Time elapsed: {time.time() - t0:.2f}s")
    grafo.get_graph_stats()

    # 3. Topological cleaning
    print(f"\n[3] TOPOLOGICAL CLEANING")

    #Tips
    print("> Processing Tips...")
    tips_removed = grafo.process_tips()
    print(f"Tips removed: {tips_removed}")

    #Bubbles
    print("> Processing Bubbles...")
    bubbles_removed = grafo.process_bubbles(cov_threshold=4)
    print(f"Bubbles resolved: {bubbles_removed}")

    #Islands
    print("> Processing Isolated Islands...")
    grafo.process_short_islands()

    print("\n[CLEAN GRAPH STATUS]")
    grafo.get_graph_stats()

    # 4. Assembly
    print(f"\n[4] ASSEMBLY")

    t0 = time.time()
    contigs_tuples = grafo.assemble_contigs(fork_in_hub=True)

    print(f"> Time elapsed: {time.time() - t0:.2f}s")
    print(f"> Raw contigs generated: {len(contigs_tuples)}")

    # 5. Final output and Statistics
    print(f"\n[5] WRITING OUTPUT")

    try:
        saved_lengths = grafo.write_contigs(OUTPUT_FILE, verbose = True)

        if saved_lengths:
            num_contigs = len(saved_lengths)
            max_len = max(saved_lengths)
            n50_val = grafo.n50(saved_lengths)

            print(f"\nSUCCESS! Assembly saved to: {OUTPUT_FILE}")
            print("-" * 40)
            print(f"Total Contigs: {num_contigs}")
            print(f"Max Length:    {max_len} bp")
            print(f"N50:           {n50_val}")
            print("-" * 40)

        else:
            print("WARNING: No contigs saved. Check length filters or assembly logic.")

    except Exception as e:
        print(f"Error writing file: {e}")

    print(f"\nTotal Pipeline Duration: {time.time() - start_global:.2f}s")


if __name__ == "__main__":
    main()