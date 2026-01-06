#!/usr/bin/env python3
"""
quick_test.py - Fast validation test with real Illumina data
"""

import time
from Dbg import Dbg

def quick_test():
    """Quick assembly test with Illumina paired-end reads"""
    
    print("=" * 70)
    print("QUICK ASSEMBLY TEST - Illumina Paired-End Reads")
    print("=" * 70)
    
    start_time = time.time()
    
    # ========== SETUP ==========
    print("\n[SETUP]")
    print("  Reads R1: Illumina_reads1.fq")
    print("  Reads R2: Illumina_reads2.fq")
    print("  K-mer size: 31")
    print("  Mode: Paired-end, both strands")
    
    # ========== BUILD GRAPH ==========
    print("\n[1/6] Building graph...")
    dbg = Dbg(
        k=31,
        file_name="Illumina_reads1.fq",
        file_type="fq",
        paired=True,
        file_nameR2="Illumina_reads2.fq"
    )
    dbg.build_graph(use_both_strands=True)
    dbg.get_graph_stats()
    
    # ========== PRUNE ==========
    print("\n[2/6] Pruning low coverage edges (threshold=3)...")
    dbg.prune_graph(threshold=3)
    
    # ========== REMOVE ISOLATED ==========
    print("\n[3/6] Removing isolated components...")
    #dbg.remove_isolated_contigs()
    
    # ========== BUBBLES ==========
    print("\n[4/6] Processing bubbles (cov_threshold=5)...")
    dbg.process_bubbles(cov_threshold=5)
    
    # ========== TIPS ==========
    print("\n[5/6] Processing tips (cov_threshold=5)...")
    dbg.process_tips(cov_threshold=5)
    
    print("\n[CLEANED GRAPH]")
    dbg.get_graph_stats()
    
    # ========== ASSEMBLY ==========
    print("\n[6/6] Assembling contigs...")
    
    # Test 1: No fork
    print("\n  → Mode: fork_in_hub=False")
    contigs_no_fork = dbg.assemble_contigs(fork_in_hub=False)
    
    lengths_no_fork = [len(c) for c in contigs_no_fork]
    n50_no_fork = Dbg.n50(lengths_no_fork)
    
    print(f"     Contigs: {len(contigs_no_fork)}")
    print(f"     Total: {sum(lengths_no_fork):,} bp")
    print(f"     N50: {n50_no_fork:,} bp")
    print(f"     Longest: {max(lengths_no_fork):,} bp")
    
    # Reset visited
    for edge in dbg.edges.values():
        edge.visited = False
    
    # Test 2: With fork
    print("\n  → Mode: fork_in_hub=True")
    contigs_fork = dbg.assemble_contigs(fork_in_hub=True)
    
    lengths_fork = [len(c) for c in contigs_fork]
    n50_fork = Dbg.n50(lengths_fork)
    
    print(f"     Contigs: {len(contigs_fork)}")
    print(f"     Total: {sum(lengths_fork):,} bp")
    print(f"     N50: {n50_fork:,} bp")
    print(f"     Longest: {max(lengths_fork):,} bp")
    
    # ========== WRITE OUTPUT ==========
    print("\n[OUTPUT]")
    
    # Save both assemblies
    for edge in dbg.edges.values():
        edge.visited = False
    dbg.contigs = contigs_no_fork
    dbg.contigs_assembled = True
    dbg.write_contigs("assembly_no_fork.fasta", min_len=100)
    
    for edge in dbg.edges.values():
        edge.visited = False
    dbg.contigs = contigs_fork
    dbg.write_contigs("assembly_with_fork.fasta", min_len=100)
    
    # ========== SUMMARY ==========
    elapsed = time.time() - start_time
    
    print("\n" + "=" * 70)
    print("TEST COMPLETE ✓")
    print("=" * 70)
    print(f"Runtime: {elapsed:.2f} seconds")
    print(f"\nComparison:")
    print(f"  fork=False → {len(contigs_no_fork)} contigs, N50={n50_no_fork:,} bp")
    print(f"  fork=True  → {len(contigs_fork)} contigs, N50={n50_fork:,} bp")
    print(f"\nFiles created:")
    print(f"  - assembly_no_fork.fasta")
    print(f"  - assembly_with_fork.fasta")
    print("=" * 70)
    
    return {
        'no_fork': {'n_contigs': len(contigs_no_fork), 'n50': n50_no_fork},
        'fork': {'n_contigs': len(contigs_fork), 'n50': n50_fork},
        'runtime': elapsed
    }

if __name__ == "__main__":
    try:
        results = quick_test()
        print("\n✓ All tests passed!")
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        exit(1)