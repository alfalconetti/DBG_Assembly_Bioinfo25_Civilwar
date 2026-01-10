from dbg import Dbg
def run_assembly_pipeline(
    k, 
    file_name, 
    file_type, 
    # Parametri Input
    paired=False, 
    file_nameR2=None, 
    reverse_r2=False,
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
    Esegue l'intera pipeline di de Bruijn Graph assembly.
    """
    
    # 1. Inizializzazione della classe Dbg
    if verbose: 
        print(f"--- DBG Initialization (k={k}) ---")
    dbg = Dbg(k, file_name, file_type, paired, file_nameR2, reverse_r2)

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

# Esempio di utilizzo su dati paired:
if __name__ == "__main__":
    run_assembly_pipeline(k=31, file_name="test_data/Illumina_reads1.fq", paired = True, file_type="fastq", quality_check=True, file_nameR2 = "test_data/Illumina_reads2.fq")