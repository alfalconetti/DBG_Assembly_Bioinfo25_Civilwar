from src.pipeline_function import run_assembly_pipeline

# Esempio di utilizzo su dati paired:
if __name__ == "__main__":
    import os
    test_r1 = "test_data/Illumina_reads1.fq"
    test_r2 = "test_data/Illumina_reads2.fq"
    if os.path.exists(test_r1) and os.path.exists(test_r2):
        print("Running test assembly")
        run_assembly_pipeline(k=31, file_name=test_r1, paired = True, file_type="fastq", quality_check=True, file_nameR2 = test_r2)

    else:
        print("Test data not found.")