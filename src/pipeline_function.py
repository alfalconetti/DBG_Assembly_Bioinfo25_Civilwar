from src.dbg import Dbg

def run_assembly_pipeline(
    k, 
    file_name, 
    file_type, 
    # Parametri Input
    paired=False, 
    file_nameR2=None, 
    reverse_r2=False,
    offset = 33,
    # Parametri Costruzione e Qualità
    quality_check=False, 
    quality_threshold=30, 
    use_both_strands=True,
    # Parametri Cleaning
    prune_threshold=3, # Più è basso più conservo
    bubble_cov_threshold=4, # Più è basso più conservo bolle
    bubble_ratio_threshold=0.1,
    tip_cov_threshold=4, # Più è basso più conservo tips
    remove_islands=True,
    # Parametri Output
    output_fasta="assembled_contigs.fasta",
    min_contig_len=0,
    fork_in_hub=False,
    verbose=True
):
    """
    Executes the complete De Bruijn Graph (DBG) assembly pipeline.

    This function orchestrates the entire workflow, from raw read parsing to 
    the generation of topologically classified contigs in FASTA format. It 
    handles graph construction, graph pruning, and structural variant 
    resolution (tips and bubbles).

    :param k: K-mer size for graph construction (must be an odd integer).
    :type k: int
    :param file_name: Path to the primary sequencing file (R1 if paired).
    :type file_name: str
    :param file_type: Format of the input file (e.g., 'fastq', 'fasta').
    :type file_type: str
    :param paired: Whether the input consists of paired-end reads, defaults to False.
    :type paired: bool
    :param file_nameR2: Path to the secondary sequencing file (R2), defaults to None.
    :type file_nameR2: str or None
    :param reverse_r2: Whether R2 reads are already reverse complemented, defaults to False.
    :type reverse_r2: bool
    :param offset:The Phred quality offset used by the sequencing platform. Defaults to 33. Can be 33 or 64
    :type offset: int
    :param quality_check: Enable filtering k-mers by base quality, defaults to False.
    :type quality_check: bool
    :param quality_threshold: Minimum Phred score to retain a k-mer, defaults to 30.
    :type quality_threshold: int
    :param use_both_strands: Whether to include reverse complements in the graph, defaults to True.
    :type use_both_strands: bool
    :param prune_threshold: Frequency threshold to remove low-coverage k-mers, defaults to 3.
    :type prune_threshold: int
    :param bubble_cov_threshold: Max coverage for a branch to be resolved as a bubble, defaults to 4.
    :type bubble_cov_threshold: int
    :param bubble_ratio_threshold: Max coverage ratio between bubble branches, defaults to 0.1.
    :type bubble_ratio_threshold: float
    :param tip_cov_threshold: Max coverage for a dead-end branch to be removed, defaults to 4.
    :type tip_cov_threshold: int
    :param remove_islands: Whether to discard small disconnected components, defaults to True.
    :type remove_islands: bool
    :param output_fasta: Name of the resulting FASTA file, defaults to "assembled_contigs.fasta".
    :type output_fasta: str
    :param min_contig_len: Minimum length of contigs to be saved, defaults to 0.
    :type min_contig_len: int
    :param fork_in_hub: If True, walking continues through branching hubs up to one level of depth, defaults to False.
    :type fork_in_hub: bool
    :param verbose: If True, prints progress and statistics to stdout, defaults to True.
    :type verbose: bool

    :return: A tuple containing (Dbg object, list of contig lengths, N50 value).
    :rtype: (Dbg, list[int], int)
    """
    
    # 1. Inizializzazione della classe Dbg
    if verbose: 
        print(f"--- DBG Initialization (k={k}) ---")
    dbg = Dbg(k, file_name, file_type, paired, file_nameR2, reverse_r2, offset)

    # 2. Costruzione del Grafo
    if verbose: 
        print("--- Buildining Graph... ---")
    
    dbg.build_graph(
        quality_check=quality_check, 
        quality_threshold=quality_threshold, 
        use_both_strands=use_both_strands,
        verbose = verbose
    )
    
    # 3. Pruning iniziale (rimozione k-meri a bassa frequenza/errori)
    if verbose:
        print(f"--- Pruning graph (threshold={prune_threshold}) ---")
    dbg.prune_graph(threshold=prune_threshold, verbose = verbose)

    # 4. Pulizia Topologica: Tips (Dead ends)
    if verbose: 
        print("--- Processing Tips ---")
    processed_tips = dbg.process_tips(cov_threshold=tip_cov_threshold, verbose = verbose)
    if verbose:
        print(f"Removed {processed_tips} tips")

    # 5. Pulizia Topologica: Bubbles (Varianti/Errori interni)
    if verbose: 
        print("--- Processing Bubbles ---")
    processed_bubbles = dbg.process_bubbles(cov_threshold=bubble_cov_threshold, ratio_threshold=bubble_ratio_threshold, verbose = verbose)
    if verbose:
        print(f"Removed {processed_bubbles} bubbles")

    # 6. Rimozione Isole corte (opzionale)
    if remove_islands:
        if verbose: 
            print("--- Removing short islands ---")
        dbg.process_short_islands(verbose = verbose)
        

    # 7. Statistiche pre-assembly
    if verbose:
        print("\n--- Graph final stats before assembly ---")
        dbg.get_graph_stats()

    # 8. Assembly dei Contigs
    if verbose: 
        print(f"\n--- Assembling contigs (fork_in_hub={fork_in_hub}) ---")
    dbg.assemble_contigs(fork_in_hub=fork_in_hub)

    # 9. Scrittura Output e Calcolo Metriche
    if verbose: 
        print(f"--- Writing output file: {output_fasta} ---")
    lengths = dbg.write_contigs(
        output_filename=output_fasta, 
        min_len=min_contig_len, 
        verbose=verbose
    )

    # 10. Risultati Finali
    if lengths:
        n50_val = dbg.n50(lengths)
        total_bp = sum(lengths)
        if verbose:
            print("\n" + "="*30)
            print(f"ASSEMBLY COMPLETED")
            print(f"Generated contigs: {len(lengths)}")
            print(f"Total bases: {total_bp} bp")
            print(f"N50: {n50_val} bp")
            print("="*30)
        return dbg, lengths, n50_val
    else:
        print("WARNING: no contig generated")
        return dbg, [], 0
