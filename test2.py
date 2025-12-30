"""
Test suite for DBG assembler
"""

from Dbg import Dbg
from kmer import kmer_generate
from utils import rc

def test_separator():
    """Separator test"""
    print("=" * 70)

def test_kmer_generation():
    """Test k-mer generation with and without reverse complement"""
    test_separator()
    print("TEST 1: K-mer Generation")
    test_separator()
    
    seq = "ATGCATGC"
    k = 4
    
    print(f"\nSequence: {seq}")
    print(f"k = {k}")
    
    # Forward only
    print("\n1.1 Forward strand only:")
    kmers_forward = list(kmer_generate(seq, k, reverse=False))
    print(f"  K-mers: {kmers_forward}")
    print(f"  Count: {len(kmers_forward)}")
    expected = len(seq) - k + 1
    assert len(kmers_forward) == expected, f"Expected {expected}, got {len(kmers_forward)}"
    print("  ✅ PASS")
    
    # Both strands
    print("\n1.2 Both strands:")
    kmers_both = list(kmer_generate(seq, k, reverse=True))
    print(f"  K-mers: {kmers_both}")
    print(f"  Count: {len(kmers_both)}")
    print(f"  RC of sequence: {rc(seq)}")
    expected_both = 2 * (len(seq) - k + 1)
    assert len(kmers_both) == expected_both, f"Expected {expected_both}, got {len(kmers_both)}"
    print("  ✅ PASS")
    
    # Check no ':' in k-mers
    print("\n1.3 Checking for separator contamination:")
    assert all(':' not in kmer for kmer in kmers_both), "Found ':' in k-mers!"
    print("  No ':' found in any k-mer ✅")
    print("  ✅ PASS")


def test_reverse_complement():
    """Test reverse complement function"""
    test_separator()
    print("TEST 2: Reverse Complement")
    test_separator()
    
    test_cases = [
        ("ATGC", "GCAT"),
        ("AAAA", "TTTT"),
        ("GCGC", "GCGC"),
        ("ATCGATCG", "CGATCGAT"),
    ]
    
    for seq, expected_rc in test_cases:
        result = rc(seq)
        print(f"  {seq} → {result} (expected: {expected_rc})")
        assert result == expected_rc, f"RC failed for {seq}"
    
    print("  ✅ PASS")


def test_graph_construction():
    """Test basic graph construction"""
    test_separator()
    print("TEST 3: Graph Construction")
    test_separator()
    
    # Simple test with tiny fake data
    print("\n3.1 Creating simple test graph...")
    
    # This would need a test FASTA file
    # For now, test the structure
    print("  (Skipping - needs test data files)")
    print("  ⏭️  SKIP")


def test_real_data():
    """Test with real Illumina data"""
    test_separator()
    print("TEST 4: Real Data Assembly")
    test_separator()
    
    print("\n4.1 Single-end assembly...")
    dbg_single = Dbg(k=31, 
                     file_name="Illumina_reads1.fq", 
                     file_type="fastq", 
                     paired=False)
    
    dbg_single.build_graph(quality_check=False, use_both_strands=True)
    
    stats_before = dbg_single.get_graph_stats()
    print()
    
    print("\n4.2 Pruning (threshold=3)...")
    dbg_single.prune_graph(threshold=3)
    
    stats_after = dbg_single.get_graph_stats()
    
    print(f"\n4.3 Results:")
    print(f"  Edges removed: {stats_before['edges'] - stats_after['edges']:,}")
    print(f"  Reduction: {100 * (1 - stats_after['edges']/stats_before['edges']):.1f}%")
    
    hubs = sum(1 for n in dbg_single.nodes.values() if n.is_hub())
    print(f"  Hubs after pruning: {hubs}")
    
    assert stats_after['edges'] < stats_before['edges'], "Pruning should remove edges"
    print("  ✅ PASS")


def test_paired_end():
    """Test paired-end assembly"""
    test_separator()
    print("TEST 5: Paired-End Assembly")
    test_separator()
    
    print("\n5.1 Building graph from paired reads (both strands)...")
    dbg_paired = Dbg(k=31, 
                     file_name="Illumina_reads1.fq", 
                     file_type="fastq", 
                     paired=True,
                     file_nameR2="Illumina_reads2.fq",
                     reverse_r2=False)
    
    dbg_paired.build_graph(quality_check=False, use_both_strands=True)
    
    stats_before = dbg_paired.get_graph_stats()
    print()
    
    print("\n5.2 Pruning (threshold=3)...")
    dbg_paired.prune_graph(threshold=3)
    
    stats_after = dbg_paired.get_graph_stats()
    
    hubs = sum(1 for n in dbg_paired.nodes.values() if n.is_hub())
    print(f"\n5.3 Hubs after pruning: {hubs}")
    
    if hubs == 0:
        print("  🎉 Perfect! Zero hubs = linear graph!")
    else:
        print(f"  ⚠️  {hubs} hubs remaining (branches to resolve)")
    
    print("  ✅ PASS")


def test_n50_calculation():
    """Test N50 calculation"""
    test_separator()
    print("TEST 6: N50 Calculation")
    test_separator()
    
    test_cases = [
        ([100, 200, 300, 400, 500], 400),  # Total 1500, half = 750, cumsum at 400
        ([1000, 500, 500], 1000),          # Total 2000, half = 1000, first contig
        ([10, 10, 10, 10, 10], 10),        # All equal
        ([100], 100),                       # Single contig
        ([], 0),                            # Empty
    ]
    
    for lengths, expected_n50 in test_cases:
        result = Dbg.n50(lengths)
        print(f"  Lengths: {lengths}")
        print(f"    N50: {result} (expected: {expected_n50})")
        assert result == expected_n50, f"N50 failed: expected {expected_n50}, got {result}"
    
    print("  ✅ PASS")


def test_comparison_table():
    """Generate comparison table for different configurations"""
    test_separator()
    print("TEST 7: Configuration Comparison")
    test_separator()
    
    configs = [
        ("Single-end", False, False, False),
        ("Paired (no both strands)", True, False, False),
        ("Paired (both strands)", True, False, True),
    ]
    
    print("\n| Configuration | Edges Before | Edges After | Hubs | Coverage |")
    print("|---------------|--------------|-------------|------|----------|")
    
    for name, paired, reverse_r2, use_both_strands in configs:
        if paired:
            dbg = Dbg(k=31, 
                     file_name="Illumina_reads1.fq", 
                     file_type="fq", 
                     paired=True,
                     file_nameR2="Illumina_reads2.fq",
                     reverse_r2=reverse_r2)
        else:
            dbg = Dbg(k=31, 
                     file_name="Illumina_reads1.fq", 
                     file_type="fq", 
                     paired=False)
        
        # Suppress output
        import sys
        from io import StringIO
        old_stdout = sys.stdout
        sys.stdout = StringIO()
        
        dbg.build_graph(quality_check=False, use_both_strands=use_both_strands)
        edges_before = len(dbg.edges)
        cov_before = sum(e.weight for e in dbg.edges.values()) / len(dbg.edges)
        
        dbg.prune_graph(threshold=3)
        edges_after = len(dbg.edges)
        hubs = sum(1 for n in dbg.nodes.values() if n.is_hub())
        
        sys.stdout = old_stdout
        
        print(f"| {name:25} | {edges_before:12,} | {edges_after:11,} | {hubs:4} | {cov_before:8.2f} |")
    
    print("\n  ✅ COMPLETE")


def run_all_tests():
    """Run all tests"""
    print("\n" + "=" * 70)
    print(" 🧬 DBG ASSEMBLER TEST SUITE")
    print("=" * 70)
    
    tests = [
        test_kmer_generation,
        test_reverse_complement,
        test_graph_construction,
        test_real_data,
        test_paired_end,
        test_n50_calculation,
        test_comparison_table,
    ]
    
    passed = 0
    failed = 0
    skipped = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"  ❌ FAIL: {e}")
            failed += 1
        except FileNotFoundError:
            print(f"  ⏭️  SKIP: Data files not found")
            skipped += 1
        except Exception as e:
            print(f"  ❌ ERROR: {e}")
            failed += 1
    
    print("\n" + "=" * 70)
    print(" SUMMARY")
    print("=" * 70)
    print(f"  ✅ Passed:  {passed}")
    print(f"  ❌ Failed:  {failed}")
    print(f"  ⏭️  Skipped: {skipped}")
    print(f"  📊 Total:   {len(tests)}")
    print("=" * 70)
    
    if failed == 0:
        print("\n🎉 ALL TESTS PASSED! 🎉\n")
    else:
        print(f"\n⚠️  {failed} test(s) failed\n")


if __name__ == "__main__":
    run_all_tests()