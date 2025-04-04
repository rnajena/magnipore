class TestGenomic:
    # Successfully reads a SAF file with all required columns
    def test_read_saf_with_required_columns(self):
        import pandas as pd
        import tempfile
        import os
        from src.genomic import read_annot
    
        # Create a temporary SAF file with required columns
        with tempfile.NamedTemporaryFile(suffix='.saf', mode='w', delete=False) as temp_file:
            temp_file.write("GeneID\tChr\tStart\tEnd\tStrand\n")
            temp_file.write("gene1\tchr1\t100\t200\t+\n")
            temp_file.write("gene2\tchr2\t300\t400\t-\n")
            temp_file_path = temp_file.name
    
        try:
            # Read the SAF file
            result_df = read_annot(temp_file_path)
        
            # Check if the DataFrame has the expected shape and content
            assert isinstance(result_df, pd.DataFrame)
            assert result_df.shape == (2, 5)
            assert list(result_df.columns) == ["GeneID", "Chr", "Start", "End", "Strand"]
            assert result_df.iloc[0]["GeneID"] == "gene1"
            assert result_df.iloc[1]["Chr"] == "chr2"
            assert result_df.iloc[0]["Start"] == 100
            assert result_df.iloc[1]["End"] == 400
            assert result_df.iloc[1]["Strand"] == "-"
        finally:
            # Clean up the temporary file
            os.unlink(temp_file_path)

    # Handles SAF file missing required columns
    def test_read_saf_missing_required_columns(self):
        import tempfile
        import os
        import pytest
        from src.genomic import read_annot
    
        # Create a temporary SAF file with missing required columns
        with tempfile.NamedTemporaryFile(suffix='.saf', mode='w', delete=False) as temp_file:
            temp_file.write("GeneID\tChr\tStart\tEnd\n")  # Missing Strand column
            temp_file.write("gene1\tchr1\t100\t200\n")
            temp_file.write("gene2\tchr2\t300\t400\n")
            temp_file_path = temp_file.name
    
        try:
            # Attempt to read the SAF file with missing columns should raise ValueError
            with pytest.raises(ValueError) as excinfo:
                read_annot(temp_file_path)
        
            # Check that the error message mentions the missing columns
            assert "SAF file must contain columns" in str(excinfo.value)
            assert "Strand" in str(excinfo.value)
        finally:
            # Clean up the temporary file
            os.unlink(temp_file_path)

    # Successfully reads a GFF file and extracts gene IDs
    def test_read_gff_extracts_gene_ids(self):
        import pandas as pd
        import tempfile
        import os
        from src.genomic import read_annot

        # Create a temporary GFF file with sample data
        with tempfile.NamedTemporaryFile(suffix='.gff', mode='w', delete=False) as temp_file:
            temp_file.write("chr1\tsource\tgene\t100\t200\t.\t+\t.\tID=gene1;\n")
            temp_file.write("chr2\tsource\tgene\t300\t400\t.\t-\t.\tID=gene2;\n")
            temp_file_path = temp_file.name

        try:
            # Read the GFF file
            result_df = read_annot(temp_file_path)
    
            # Check if the DataFrame has the expected shape and content
            assert isinstance(result_df, pd.DataFrame)
            assert result_df.shape == (2, 5)
            assert list(result_df.columns) == ["GeneID", "Chr", "Start", "End", "Strand"]
            assert result_df.iloc[0]["GeneID"] == "gene1"
            assert result_df.iloc[1]["Chr"] == "chr2"
            assert result_df.iloc[0]["Start"] == 100
            assert result_df.iloc[1]["End"] == 400
            assert result_df.iloc[1]["Strand"] == "-"
        finally:
            # Clean up the temporary file
            os.unlink(temp_file_path)

    # Successfully reads a GTF file and extracts gene IDs
    def test_read_gtf_and_extract_gene_ids(self):
        import pandas as pd
        import tempfile
        import os
        from src.genomic import read_annot

        # Create a temporary GTF file with sample data
        with tempfile.NamedTemporaryFile(suffix='.gtf', mode='w', delete=False) as temp_file:
            temp_file.write("chr1\tsource\tgene\t100\t200\t.\t+\t.\tgene_id \"gene1\";\n")
            temp_file.write("chr2\tsource\tgene\t300\t400\t.\t-\t.\tgene_id \"gene2\";\n")
            temp_file_path = temp_file.name

        try:
            # Read the GTF file
            result_df = read_annot(temp_file_path)
    
            # Check if the DataFrame has the expected shape and content
            assert isinstance(result_df, pd.DataFrame)
            assert result_df.shape == (2, 5)
            assert list(result_df.columns) == ["GeneID", "Chr", "Start", "End", "Strand"]
            assert result_df.iloc[0]["GeneID"] == "gene1"
            assert result_df.iloc[1]["Chr"] == "chr2"
            assert result_df.iloc[0]["Start"] == 100
            assert result_df.iloc[1]["End"] == 400
            assert result_df.iloc[1]["Strand"] == "-"
        finally:
            # Clean up the temporary file
            os.unlink(temp_file_path)

    # Function returns True when genomic_row and magni_row match on Chr/ref and Strand
    def test_returns_true_when_rows_match(self):
        import pandas as pd
        from src.genomic import f_match
    
        # Create test data
        genomic_row = pd.Series({'Chr': 'chr1', 'Strand': '+'})
        magni_row = pd.Series({'ref_1': 'chr1', 'strand': '+', 'ref_2': 'chr2'})
    
        # Test with sample 1
        result = f_match(genomic_row, magni_row, '1')
    
        # Assert
        assert result is True

    # Successfully reads a valid TSV file and returns a pandas DataFrame
    def test_read_valid_tsv_file(self):
        import pandas as pd
        import tempfile
        import os
        from src.genomic import read_magnipore
    
        # Create a temporary TSV file
        with tempfile.NamedTemporaryFile(suffix='.tsv', delete=False) as tmp:
            tmp.write(b"col1\tcol2\tcol3\n1\t2\t3\n4\t5\t6\n")
            tmp_path = tmp.name
    
        try:
            # Call the function
            result = read_magnipore(tmp_path)
        
            # Verify the result
            assert isinstance(result, pd.DataFrame)
            assert result.shape == (2, 3)
            assert list(result.columns) == ["col1", "col2", "col3"]
            assert result.iloc[0, 0] == 1
            assert result.iloc[1, 2] == 6
        finally:
            # Clean up
            os.unlink(tmp_path)

    # Handles files with missing values
    def test_read_tsv_with_missing_values(self):
        import pandas as pd
        import tempfile
        import os
        from src.genomic import read_magnipore
    
        # Create a temporary TSV file with missing values
        with tempfile.NamedTemporaryFile(suffix='.tsv', delete=False) as tmp:
            tmp.write(b"col1\tcol2\tcol3\n1\t\t3\n4\t5\t\n")
            tmp_path = tmp.name
    
        try:
            # Call the function
            result = read_magnipore(tmp_path)
        
            # Verify the result
            assert isinstance(result, pd.DataFrame)
            assert result.shape == (2, 3)
            assert pd.isna(result.iloc[0, 1])
            assert pd.isna(result.iloc[1, 2])
            assert result.iloc[0, 2] == 3
            assert result.iloc[1, 0] == 4
        finally:
            # Clean up
            os.unlink(tmp_path)

    # Correctly identifies genes that overlap with magnipore positions
    def test_identifies_overlapping_genes(self):
        from src.genomic import get_genes
        import pandas as pd
        # Prepare test data
        annot = pd.DataFrame({
            "Chr": ["chr1", "chr1", "chr2"],
            "Strand": ["+", "-", "+"],
            "GeneID": ["gene1", "gene2", "gene3"],
            "Start": [100, 500, 200],
            "End": [300, 700, 400]
        })
    
        magnipore = pd.DataFrame({
            "ref_1": ["chr1", "chr1", "chr2"],
            "strand": ["+", "-", "+"],
            "pos_1": [200, 600, 300],
            "base_1": ["A", "C", "G"],
            "motif_1": ["AAAT", "CCCT", "GGGT"]
        })
    
        result = get_genes(annot, magnipore, "1")
    
        # Assert results
        assert len(result) == 3
        assert result.iloc[0]["GeneID"] == "gene1"
        assert result.iloc[0]["Magnipore"] == 200
        assert result.iloc[0]["Geneposition"] == 100
        assert result.iloc[1]["GeneID"] == "gene2"
        assert result.iloc[1]["Magnipore"] == 600
        assert result.iloc[2]["GeneID"] == "gene3"
        assert result.iloc[2]["Magnipore"] == 300

    # Parsing valid arguments with all required parameters (annot, magnipore, outfile)
    def test_parse_with_valid_required_arguments(self):
        import sys
        from unittest.mock import patch
        from src.genomic import parse
    
        test_args = ['magnipore', 'test.saf', 'test.magnipore', 'output.tsv']
        with patch.object(sys, 'argv', test_args):
            args = parse()
        
        assert args.annot == 'test.saf'
        assert args.magnipore == 'test.magnipore'
        assert args.outfile == 'output.tsv'
        assert args.sample == 1  # Default value

    # Missing required arguments (annot, magnipore, outfile)
    def test_parse_with_missing_required_arguments(self):
        import sys
        import pytest
        from unittest.mock import patch
        from src.genomic import parse
    
        # Test with missing arguments
        test_args = ['magnipore']
        with patch.object(sys, 'argv', test_args):
            with pytest.raises(SystemExit):
                parse()
            
        # Test with partial arguments
        test_args = ['magnipore', 'test.saf']
        with patch.object(sys, 'argv', test_args):
            with pytest.raises(SystemExit):
                parse()