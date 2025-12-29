"""
COMPREHENSIVE DBG ASSEMBLY TEST
Test completo del De Bruijn Graph Assembler
"""

from Dbg import Dbg

def print_separator(title):
    """Print a nice separator"""
    print("\n" + "=" * 80)
    print(f" {title}")
    print("=" * 80)


def test_graph_stats(dbg):
    """Calculate and print detailed graph statistics"""
    # Basic stats
    num_nodes = len(dbg.nodes)
    num_edges = len(dbg.edges)
    
    # Coverage stats
    weights = [edge.weight for edge in dbg.edges.values()]
    total_coverage = sum(weights)
    avg_coverage = total_coverage / num_edges if num_edges > 0 else 0
    min_weight = min(weights) if weights else 0
    max_weight = max(weights) if weights else 0
    median_weight = sorted(weights)[len(weights)//2] if weights else 0
    
    # Node degree analysis
    isolated = sum(1 for n in dbg.nodes.values() 
                   if len(n.incoming) == 0 and len(n.outgoing) == 0)
    start_nodes = sum(1 for n in dbg.nodes.values() 
                      if len(n.incoming) == 0 and len(n.outgoing) > 0)
    end_nodes = sum(1 for n in dbg.nodes.values() 
                    if len(n.incoming) > 0 and len(n.outgoing) == 0)
    hubs = sum(1 for n in dbg.nodes.values() if n.is_hub())
    
    # Weight distribution
    from collections import Counter
    weight_dist = Counter(weights)
    
    print(f"\n📊 GRAPH STATISTICS:")
    print(f"  Nodes:              {num_nodes:,}")
    print(f"  Edges:              {num_edges:,}")
    print(f"  Total coverage:     {total_coverage:,}")
    print(f"  Average coverage:   {avg_coverage:.2f}x")
    
    print(f"\n📈 EDGE WEIGHT DISTRIBUTION:")
    print(f"  Min weight:         {min_weight}")
    print(f"  Max weight:         {max_weight}")
    print(f"  Median weight:      {median_weight}")
    print(f"  Weight = 1:         {weight_dist.get(1, 0):,} edges")
    print(f"  Weight = 2:         {weight_dist.get(2, 0):,} edges")
    print(f"  Weight = 3:         {weight_dist.get(3, 0):,} edges")
    print(f"  Weight >= 4:        {sum(c for w, c in weight_dist.items() if w >= 4):,} edges")
    
    print(f"\n🔗 NODE ANALYSIS:")
    print(f"  Isolated nodes:     {isolated:,}")
    print(f"  Start nodes:        {start_nodes:,}")
    print(f"  End nodes:          {end_nodes:,}")
    print(f"  Hub nodes:          {hubs:,}")
    
    return {
        'nodes': num_nodes,
        'edges': num_edges,
        'coverage': avg_coverage,
        'hubs': hubs,
        'isolated': isolated
    }


def comprehensive_test():
    """Run comprehensive test of the DBG assembler"""
    
    print_separator("🧬 COMPREHENSIVE DBG ASSEMBLY TEST")
    
    # ==================== TEST 1: Single-end ====================
    print_separator("TEST 1: Single-End Assembly")
    
    print("\n1.1 Building graph from R1 only...")
    dbg_single = Dbg(k=31, 
                     file_name="Illumina_reads1.fq", 
                     file_type="fq", 
                     paired=False)
    
    dbg_single.build_graph(quality_check=False)
    stats_single_before = test_graph_stats(dbg_single)
    
    print("\n1.2 Pruning graph (threshold=2)...")
    dbg_single.prune_graph(threshold=2)
    stats_single_after = test_graph_stats(dbg_single)
    
    print(f"\n1.3 Pruning results:")
    print(f"  Edges removed:      {stats_single_before['edges'] - stats_single_after['edges']:,}")
    print(f"  Reduction:          {100 * (1 - stats_single_after['edges']/stats_single_before['edges']):.1f}%")
    
    
    # ==================== TEST 2: Paired-end WITHOUT RC ====================
    print_separator("TEST 2: Paired-End Assembly (NO Reverse Complement)")
    
    print("\n2.1 Building graph from R1 + R2 (no RC)...")
    dbg_paired_no_rc = Dbg(k=31, 
                           file_name="Illumina_reads1.fq", 
                           file_type="fq", 
                           paired=True,
                           file_nameR2="Illumina_reads2.fq",
                           reverse_r2=False)
    
    dbg_paired_no_rc.build_graph(quality_check=False)
    stats_paired_no_rc_before = test_graph_stats(dbg_paired_no_rc)
    
    print("\n2.2 Pruning graph (threshold=2)...")
    dbg_paired_no_rc.prune_graph(threshold=2)
    stats_paired_no_rc_after = test_graph_stats(dbg_paired_no_rc)
    
    print(f"\n2.3 Pruning results:")
    print(f"  Edges removed:      {stats_paired_no_rc_before['edges'] - stats_paired_no_rc_after['edges']:,}")
    print(f"  Reduction:          {100 * (1 - stats_paired_no_rc_after['edges']/stats_paired_no_rc_before['edges']):.1f}%")
    
    
    # ==================== TEST 3: Paired-end WITH RC ====================
    print_separator("TEST 3: Paired-End Assembly (WITH Reverse Complement)")
    
    print("\n3.1 Building graph from R1 + R2 (with RC)...")
    dbg_paired_rc = Dbg(k=31, 
                        file_name="Illumina_reads1.fq", 
                        file_type="fq", 
                        paired=True,
                        file_nameR2="Illumina_reads2.fq",
                        reverse_r2=True)
    
    dbg_paired_rc.build_graph(quality_check=False)
    stats_paired_rc_before = test_graph_stats(dbg_paired_rc)
    
    print("\n3.2 Pruning graph (threshold=2)...")
    dbg_paired_rc.prune_graph(threshold=2)
    stats_paired_rc_after = test_graph_stats(dbg_paired_rc)
    
    print(f"\n3.3 Pruning results:")
    print(f"  Edges removed:      {stats_paired_rc_before['edges'] - stats_paired_rc_after['edges']:,}")
    print(f"  Reduction:          {100 * (1 - stats_paired_rc_after['edges']/stats_paired_rc_before['edges']):.1f}%")
    
    
    # ==================== TEST 4: Different thresholds ====================
    print_separator("TEST 4: Testing Different Pruning Thresholds")
    
    thresholds = [2, 3, 4, 5]
    
    print("\n4.1 Building fresh graph...")
    dbg_threshold = Dbg(k=31, 
                        file_name="Illumina_reads1.fq", 
                        file_type="fq", 
                        paired=True,
                        file_nameR2="Illumina_reads2.fq",
                        reverse_r2=False)
    dbg_threshold.build_graph(quality_check=False)
    
    print(f"\n4.2 Testing thresholds:")
    print(f"  Initial edges: {len(dbg_threshold.edges):,}")
    
    for t in thresholds:
        # Create fresh copy for each threshold
        dbg_test = Dbg(k=31, 
                       file_name="Illumina_reads1.fq", 
                       file_type="fq", 
                       paired=True,
                       file_nameR2="Illumina_reads2.fq",
                       reverse_r2=False)
        dbg_test.build_graph(quality_check=False)
        
        edges_before = len(dbg_test.edges)
        dbg_test.prune_graph(threshold=t)
        edges_after = len(dbg_test.edges)
        hubs_after = sum(1 for n in dbg_test.nodes.values() if n.is_hub())
        
        print(f"  Threshold {t}: {edges_after:,} edges remaining "
              f"({100 * edges_after/edges_before:.1f}%), {hubs_after} hubs")
    
    
    # ==================== COMPARISON SUMMARY ====================
    print_separator("📋 COMPARISON SUMMARY")
    
    print("\n🔍 BEFORE PRUNING:")
    print(f"  Single-end:         {stats_single_before['edges']:,} edges, "
          f"{stats_single_before['hubs']} hubs")
    print(f"  Paired (no RC):     {stats_paired_no_rc_before['edges']:,} edges, "
          f"{stats_paired_no_rc_before['hubs']} hubs")
    print(f"  Paired (with RC):   {stats_paired_rc_before['edges']:,} edges, "
          f"{stats_paired_rc_before['hubs']} hubs")
    
    print("\n✨ AFTER PRUNING (threshold=2):")
    print(f"  Single-end:         {stats_single_after['edges']:,} edges, "
          f"{stats_single_after['hubs']} hubs")
    print(f"  Paired (no RC):     {stats_paired_no_rc_after['edges']:,} edges, "
          f"{stats_paired_no_rc_after['hubs']} hubs ⭐")
    print(f"  Paired (with RC):   {stats_paired_rc_after['edges']:,} edges, "
          f"{stats_paired_rc_after['hubs']} hubs")
    
    print("\n🏆 RECOMMENDATION:")
    if stats_paired_no_rc_after['hubs'] < stats_paired_rc_after['hubs']:
        print(f"  ✅ Use paired-end WITHOUT reverse complement")
        print(f"     Cleaner graph: {stats_paired_no_rc_after['hubs']} hubs vs "
              f"{stats_paired_rc_after['hubs']} hubs")
    else:
        print(f"  ✅ Use paired-end WITH reverse complement")
        print(f"     Cleaner graph: {stats_paired_rc_after['hubs']} hubs vs "
              f"{stats_paired_no_rc_after['hubs']} hubs")
    
    print("\n💡 NEXT STEPS:")
    print("  1. Confirm reverse complement strategy with professor")
    print("  2. Implement contig assembly from graph")
    print("  3. Write contigs to FASTA format")
    print("  4. Calculate assembly metrics (N50, total length)")
    
    print_separator("✅ TEST COMPLETED")
    

if __name__ == "__main__":
    try:
        comprehensive_test()
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()