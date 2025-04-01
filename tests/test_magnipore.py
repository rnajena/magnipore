import os
import numpy as np
import pandas as pd
import pytest
from Bio import SeqIO
from src import magnipore
import multiprocessing

multiprocessing.set_start_method("forkserver", force=True)  # Alternative: "spawn"

class TestMagnipore:
    @pytest.fixture(scope="module")
    def test_data(self):
        """Load test data from files and return a dictionary of necessary objects."""
        magnipore.init_Logger(None)
        this_file_dir = os.path.dirname(__file__)

        # File paths
        ref1 = os.path.join(this_file_dir, 'sample1.fa')
        red1 = os.path.join(this_file_dir, 'sample1.red')
        ref2 = os.path.join(this_file_dir, 'sample2.fa')
        red2 = os.path.join(this_file_dir, 'sample2.red')
        alignment_path = os.path.join(this_file_dir, 'alignment', 'sample1_sample2.aln')

        # Read sequences and red files
        seq_dict = {**SeqIO.to_dict(SeqIO.parse(ref1, "fasta")), **SeqIO.to_dict(SeqIO.parse(ref2, "fasta"))}
        red1DF = pd.read_csv(red1, sep='\t')
        red2DF = pd.read_csv(red2, sep='\t')

        # Process alignment
        mapping_dict, unaligned_dict, aln_dict, seq_dict = magnipore.getMapping(alignment_path, this_file_dir, 'sample1', 'sample2')
        red1_dict = magnipore.read_red_file(red1, list(seq_dict.values())[0])
        red2_dict = magnipore.read_red_file(red2, list(seq_dict.values())[-1])

        return {
            "this_file_dir": this_file_dir,
            "red1DF": red1DF,
            "red2DF": red2DF,
            "red1_dict": red1_dict,
            "red2_dict": red2_dict,
            "seq_dict": seq_dict,
            "aln_dict": aln_dict,
            "mapping_dict": mapping_dict,
            "unaligned_dict": unaligned_dict,
        }

    def test_read_red_file(self, test_data):
        """Test the `read_red_file` function by validating against reference data."""
        red1DF = test_data["red1DF"]
        red1_dict = test_data["red1_dict"]

        for i, pos in enumerate(red1_dict):
            for j, red in enumerate(pos):
                assert np.isclose(red.mean, red1DF.loc[(2*i)+j, 'signal_mean']), "Mismatch in signal_mean"
                assert np.isclose(red.std, red1DF.loc[(2*i)+j, 'signal_std']), "Mismatch in signal_std"
                assert np.isclose(red.data_density, red1DF.loc[(2*i)+j, 'data_density']), "Mismatch in data_density"
                assert np.isclose(red.n_datapoints, red1DF.loc[(2*i)+j, 'n_datapoints']), "Mismatch in n_datapoints"
                assert np.isclose(red.contained_datapoints, red1DF.loc[(2*i)+j, 'contained_datapoints']), "Mismatch in contained_datapoints"
                assert np.isclose(red.n_segments, red1DF.loc[(2*i)+j, 'n_segments']), "Mismatch in n_segments"
                assert np.isclose(red.contained_segments, red1DF.loc[(2*i)+j, 'contained_segments']), "Mismatch in contained_segments"
                assert np.isclose(red.n_reads, red1DF.loc[(2*i)+j, 'n_reads']), "Mismatch in n_reads"

    def test_mapping(self, test_data):
        """Test the `getMapping` function by validating the output against expected results."""
        expected_mapping = {
            0: (0, 0), 1: (1, 1), 2: (2, 2), 3: (3, 3), 4: (4, 4), 10: (5, 10),
            11: (6, 11), 12: (7, 12), 13: (8, 13), 14: (9, 14), 15: (10, 15),
            16: (11, 16), 17: (12, 17), 18: (13, 18), 19: (14, 19), 20: (15, 20),
            21: (16, 21)
        }
        expected_unaligned = {
            'sample1': [(5, 'T'), (6, 'T'), (7, 'T'), (8, 'T'), (9, 'T'), (22, 'T')],
            'sample2': []
        }

        assert test_data["mapping_dict"] == expected_mapping, "Mapping dictionary mismatch"
        assert test_data["unaligned_dict"] == expected_unaligned, "Unaligned positions mismatch"

        # Handles single sequence input by creating identity mapping
    def test_handles_single_sequence_with_identity_mapping(self, mocker):
        # Mock the Logger
        # mock_logger = mocker.patch('src.magnipore.LOGGER')
    
        # Create temporary alignment file
        import tempfile
        from os.path import join
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio import SeqIO
    
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test alignment file with single sequence
            seq1 = SeqRecord(Seq("ACGTACGT"), id="seq1")
        
            align_file = join(tmpdir, "test_single_seq.fasta")
            with open(align_file, "w") as f:
                SeqIO.write([seq1], f, "fasta")
        
            # Call the function
            from src.magnipore import getMapping
            mapping, unaligned, alignment, sequences = getMapping(
                alignment=align_file,
                out=tmpdir,
                l1="seq1",
                l2="none"  # This parameter is still required but not used for single sequence
            )
        
            # Verify the identity mapping
            expected_mapping = {i: (i, i) for i in range(8)}  # Length of "ACGTACGT" is 8
            assert mapping == expected_mapping
        
            # Verify unaligned positions (should be empty for single sequence)
            assert list(unaligned.keys())[0] == "seq1"
            assert unaligned["seq1"] == []
        
            # Verify sequences
            assert sequences["seq1"] == "ACGTACGT"
            assert len(sequences) == 1

    # Correctly maps positions between two aligned sequences, ensuring the mapping indices are within the correct range, and verifies unaligned positions.
    def test_maps_positions_between_aligned_sequences_corrected(self, mocker):
        # Mock the Logger
        # mock_logger = mocker.patch('src.magnipore.LOGGER')

        # Create temporary alignment file
        import tempfile
        from os.path import join
        from Bio.SeqRecord import SeqRecord
        from Bio.Seq import Seq
        from Bio import SeqIO

        with tempfile.TemporaryDirectory() as tmpdir:
            # Create test alignment file
            seq1 = SeqRecord(Seq("ACGT-ACGT"), id="seq1")
            seq2 = SeqRecord(Seq("ACGTA-CGT"), id="seq2")

            align_file = join(tmpdir, "test_alignment.fasta")
            with open(align_file, "w") as f:
                SeqIO.write([seq1, seq2], f, "fasta")

            # Call the function
            from src.magnipore import getMapping
            mapping, unaligned, alignment, sequences = getMapping(
                alignment=align_file,
                out=tmpdir,
                l1="seq1",
                l2="seq2"
            )

            # Verify the mapping
            assert mapping[0] == (0, 0)  # A-A
            assert mapping[1] == (1, 1)  # C-C
            assert mapping[2] == (2, 2)  # G-G
            assert mapping[3] == (3, 3)  # T-T
            assert 4 not in mapping      # - is unaligned
            assert mapping[5] == (5, 6)  # A-A
            assert mapping[6] == (6, 7)  # C-C
            assert mapping[7] == (7, 8)  # G-G

            # Verify unaligned positions
            assert unaligned["seq1"] == [(4, "A")]
            assert unaligned["seq2"] == [(4, "A")]

            # Verify sequences
            assert sequences["seq1"] == "ACGTACGT"
            assert sequences["seq2"] == "ACGTACGT"

    def test_magnipore(self, test_data):
        """Test the `magnipore` function by validating output file contents."""
        magnipore.magnipore(
            test_data["mapping_dict"], test_data["unaligned_dict"],
            test_data["seq_dict"], test_data["aln_dict"],
            test_data["red1_dict"], test_data["red2_dict"],
            'sample1', 'sample2', test_data["this_file_dir"], 5, 8
        )

        magn_file = os.path.join(test_data["this_file_dir"], 'magnipore/sample1_sample2/sample1_sample2.magnipore')
        magn = pd.read_csv(magn_file, sep='\t')

        assert not magn.empty, "Magnipore output file is empty"
        assert "td_score" in magn.columns, "Missing td_score column in magnipore output"
        assert "bayesian_p" in magn.columns, "Missing bayesian_p column in magnipore output"

    def test_stockholm(self, test_data):
        """Test if the marked stockholm file is correctly written with expected base modifications."""
        stk_file = os.path.join(test_data["this_file_dir"], 'magnipore/sample1_sample2/sample1_sample2_marked.stk')

        with open(stk_file, 'r') as f:
            for line in f:
                line = line.strip().split(' ')
                if line[0] == 'magnipore_marked_sample1':
                    assert line[-1].count('X') == 4, "Incorrect number of X marks for sample1"
                    assert line[-1].count('.') == 19, "Incorrect number of . marks for sample1"
                elif line[0] == 'magnipore_marked_sample2':
                    assert line[-1].count('X') == 4, "Incorrect number of X marks for sample2"
                    assert line[-1].count('.') == 13, "Incorrect number of . marks for sample2"
                    assert line[-1].count('-') == 6, "Incorrect number of - marks for sample2"

    # Verify that the Logger is initialized correctly with a valid file path
    def test_init_logger_with_valid_file_path(self, mocker):
        # Arrange
        from src.magnipore import init_Logger

        mock_logger = mocker.patch('src.magnipore.Logger')
        # mock_logger_instance = mock_logger.return_value

        test_file_path = "test_log.txt"

        # Act
        init_Logger(test_file_path)

        # Assert
        mock_logger.assert_called_once_with(test_file_path)

    # Initialize logger without a file path (None)
    def test_init_logger_without_file_path(self, mocker):
        # Arrange
        from src.magnipore import init_Logger

        mock_logger = mocker.patch('src.magnipore.Logger')
        # mock_logger_instance = mock_logger.return_value

        # Act
        init_Logger(None)

        # Assert
        mock_logger.assert_called_once_with(None)

    # Function correctly logs error message and exits with error code 3
    def test_logs_error_and_exits_with_code_3(self, mocker):
        # Arrange
        from src.magnipore import callbackErrorRed
        mock_logger = mocker.patch('src.magnipore.LOGGER')
        test_error = Exception("Test error message")
    
        # Act
        callbackErrorRed(test_error)
    
        # Assert
        mock_logger.error.assert_called_once_with(
            f'Error in multiprocessing red building: {test_error}', 3
        )

    # Correctly calls LOGGER.error with formatted error message
    def test_calls_logger_error_with_formatted_message(self, mocker):
        # Arrange
        from src.magnipore import callbackErrorComparison
        mock_logger = mocker.Mock()
        mocker.patch('src.magnipore.LOGGER', mock_logger)
        test_error = ValueError("Test error message")
    
        # Act
        callbackErrorComparison(test_error)
    
        # Assert
        mock_logger.error.assert_called_once_with(
            f'Error in multiprocessing magnipore signal comparison: {test_error}', 4
        )

    # Replaces all non-gap characters with dots in a sequence
    def test_replaces_non_gap_characters_with_dots(self):
        # Arrange
        from src.magnipore import reformat
        sequence = "ACGT-ACGT-ACGT"
        expected = [".", ".", ".", ".", "-", ".", ".", ".", ".", "-", ".", ".", ".", "."]
    
        # Act
        result = reformat(sequence)
    
        # Assert
        assert result == expected

    # Reconstructs Red objects from chunks while preserving their statistics
    def test_reconstructs_reds_preserving_statistics_with_t_4(self):
        import numpy as np
        from src.Red import Red

        # Create original Reds with known values
        original_reds = []
        for i in range(8):
            red = Red()
            # Add different values to each Red
            red.append(np.array([i*10 + j for j in range(5)], dtype=np.float32))
            red.get_mean_stdev()  # Calculate stats
            original_reds.append(red)

        # Divide into chunks (t=4)
        t = 4
        red_chunks = [
            [original_reds[i] for i in range(0, 8, 4)],  # [0, 4]
            [original_reds[i] for i in range(1, 8, 4)],  # [1, 5]
            [original_reds[i] for i in range(2, 8, 4)],  # [2, 6]
            [original_reds[i] for i in range(3, 8, 4)]   # [3, 7]
        ]

        # Reconstruct
        from src.magnipore import reconstruct_reds
        reconstructed = reconstruct_reds(red_chunks, len(original_reds), t)

        # Verify statistics are preserved
        for i, red in enumerate(original_reds):
            orig_mean, orig_std = red.get_mean_stdev()
            recon_mean, recon_std = reconstructed[i].get_mean_stdev()
    
            assert np.isclose(orig_mean, recon_mean)
            assert np.isclose(orig_std, recon_std)
            assert red.n == reconstructed[i].n, f"Sample count mismatch at index {i}"

    # Function correctly maps read IDs from a BAM file
    def test_read_id_mapping_with_pi_tag(self, mocker):
        # Mock AlignmentFile and its methods
        mock_alignment_file = mocker.patch('src.magnipore.AlignmentFile')
        mock_samfile = mock_alignment_file.return_value.__enter__.return_value
    
        # Create mock reads
        mock_read1 = mocker.MagicMock()
        mock_read1.query_name = "read1"
        mock_read1.has_tag.return_value = True
        mock_read1.get_tag.return_value = "processed_id1"
    
        mock_read2 = mocker.MagicMock()
        mock_read2.query_name = "read2"
        mock_read2.has_tag.return_value = False
    
        # Set up the mock to return our mock reads
        mock_samfile.fetch.return_value = [mock_read1, mock_read2]
    
        # Call the function
        from src.magnipore import getReadIdMap
        result = getReadIdMap("test.bam")
    
        # Verify the results
        assert result == {"read1": "processed_id1", "read2": "read2"}
        mock_alignment_file.assert_called_once_with("test.bam", "rb", check_sq=False)
        mock_samfile.fetch.assert_called_once_with(until_eof=True)

    # Ensures the updater_task processes a valid line correctly and updates RED objects with a correctly sized reds list
    def test_valid_line_processing_with_corrected_reds_size(self, mocker):
        from src.Red import Red
        from src.magnipore import updater_task
        # Setup
        queue = mocker.MagicMock()
        queue.get.side_effect = ["chrom\t10\t+\tA\tC\tread1\t100\t50\t0\t0", None]

        # Ensure reds list is large enough to handle the position index calculation
        reds = [[Red(), Red()] for _ in range(6)]  # 6 positions, 2 strands
        num_updaters = 2
        raw = "test_raw_file.pod5"
        read_id_map = {"read1": "signal1"}
        processed_counter = mocker.MagicMock()
        processed_counter.value = 0
        lock = mocker.MagicMock()

        # Mock read5_ont
        mock_r5 = mocker.MagicMock()
        mock_r5.getZNormSignal.return_value = np.ones(200)  # Mock signal
        mock_r5_read = mocker.patch('read5_ont.read', return_value=mock_r5)

        # Mock STRANDENCODER
        mocker.patch('src.magnipore.STRANDENCODER', {'+': 0, '-': 1})

        # Execute
        result = updater_task(queue, reds, num_updaters, raw, read_id_map, processed_counter, lock)

        # Assert
        assert result == reds
        mock_r5_read.assert_called_once_with(raw)
        mock_r5.getZNormSignal.assert_called_once_with("signal1")
        assert mock_r5.close.called

        # Check RED object was updated correctly
        red_obj = reds[10 // num_updaters][0]  # Corrected index calculation
        assert red_obj.n_reads == 1
        assert red_obj.n_segments == 1
        assert red_obj.n_datapoints == 50

        # Check counter was incremented
        lock.__enter__.assert_called_once()
        assert processed_counter.value == 1

    # Successfully processes queue items and updates Red objects with signal data, ensuring valid index access
    def test_processes_queue_and_updates_reds_with_valid_index(self, mocker):
        from src.Red import Red
        from src.magnipore import checker_task

        # Setup
        queue = mocker.MagicMock()
        queue.get.side_effect = [
            "chr1\t0\t+\tA\tC\tread1\t100\t50\t0\t0",
            None  # Termination signal
        ]

        # Mock Red objects
        red_mock = mocker.MagicMock(spec=Red)
        red_mock.get_mean_stdev.return_value = (0.0, 1.0)

        # Setup reds structure (2D list of Red objects)
        reds = [[red_mock, mocker.MagicMock(spec=Red)]]

        # Mock read5_ont
        mock_r5 = mocker.MagicMock()
        mock_r5.getZNormSignal.return_value = np.ones(200)  # Signal of all 1's
        mock_r5_module = mocker.patch('read5_ont.read', return_value=mock_r5)

        # Other parameters
        read_id_map = {"read1": "signal1"}
        processed_counter = mocker.MagicMock()
        processed_counter.value = 0
        lock = mocker.MagicMock()

        # Execute
        result = checker_task(
            queue=queue,
            reds=reds,
            num_updaters=1,
            raw="raw_file.pod5",
            read_id_map=read_id_map,
            cal_data_density=True,
            processed_counter=processed_counter,
            lock=lock
        )

        # Verify
        assert result is reds  # Should return the same reds object
        assert queue.get.call_count == 2
        mock_r5.getZNormSignal.assert_called_once_with("signal1")
        red_mock.get_mean_stdev.assert_called_once()
        red_mock.add_contained_segments.assert_called_once()
        red_mock.add_contained_datapoints.assert_called_once()
        red_mock.add_data_density.assert_called_once()
        assert processed_counter.value == 1
        mock_r5.close.assert_called_once()

    # Correctly updates the progress bar when the counter value changes
    def test_progress_bar_updates_with_counter_change(self, mocker):
        # Setup
        mock_tqdm = mocker.Mock()
        mock_counter = mocker.Mock()
        mock_counter.value = 0
        mock_lock = mocker.Mock()

        # Create a side effect that increases the counter value each time it's accessed
        counter_values = [0, 5, 10]
        mock_counter.value = counter_values[0]

        def side_effect():
            if len(counter_values) > 1:
                mock_counter.value = counter_values.pop(0)
            else:
                mock_counter.value = -1  # Terminate the loop
            return mock_lock

        # Mock the lock context manager
        mock_lock.__enter__ = mocker.Mock(side_effect=side_effect)
        mock_lock.__exit__ = mocker.Mock(return_value=None)

        # Execute
        from src.magnipore import progress_updater  # Import from the correct module
        progress_updater(mock_tqdm, mock_counter, mock_lock)

        # Verify
        assert mock_tqdm.update.call_count == 2
        mock_tqdm.update.assert_any_call(5)  # First update: 5 - 0 = 5
        mock_tqdm.update.assert_any_call(5)  # Second update: 10 - 5 = 5

    # Writes RED data to file with correct header and format
    def test_writes_red_data_with_correct_format(self, mocker, tmp_path):
        from src.Red import Red
        # Mock tqdm to avoid progress bar in tests
        mocker.patch('src.magnipore.tqdm')

        # Create test data using MagicMock with __str__ method
        red1 = mocker.MagicMock(spec=Red)
        red1.__str__.return_value = "10.12345678\t2.12345678\t0.12345678\t0.23456789\t100\t80\t5\t4\t10"

        red2 = mocker.MagicMock(spec=Red)
        red2.__str__.return_value = "15.12345678\t3.12345678\t0.22345678\t0.33456789\t200\t150\t10\t8\t20"

        # Create reds list structure (list of lists of Red objects)
        reds = [[red1], [red2]]

        # Create temporary file path
        red_file = tmp_path / "test_red.txt"

        # Call the function
        from src.magnipore import writeOutput
        writeOutput(str(red_file), reds)

        # Verify file content
        with open(red_file, 'r') as f:
            content = f.readlines()

        # Check header
        assert content[0] == 'strand\tposition\tsignal_mean\tsignal_std\tdata_density\texpected_model_density\tn_datapoints\tcontained_datapoints\tn_segments\tcontained_segments\tn_reads\n'

        # Check data rows
        assert content[1] == '+\t0\t10.12345678\t2.12345678\t0.12345678\t0.23456789\t100\t80\t5\t4\t10\n'
        assert content[2] == '+\t1\t15.12345678\t3.12345678\t0.22345678\t0.33456789\t200\t150\t10\t8\t20\n'

    # Test with two distinct normal distributions returns expected D statistic and p-value
    def test_distinct_distributions_return_expected_statistics(self):
        from src.magnipore import ks_test
        # Set random seed for reproducibility
        np.random.seed(42)
    
        # Define two distinct distributions
        dist1 = (0, 1)  # mean=0, std=1
        dist2 = (5, 1)  # mean=5, std=1
    
        # Call the function
        D, p = ks_test(dist1, dist2)
    
        # With these distributions and seed, we expect a high D value (close to 1)
        # and a very low p-value (indicating distributions are different)
        assert 0.5 < D <= 1.0, f"Expected high D statistic, got {D}"
        assert p < 0.05, f"Expected low p-value, got {p}"
    
        # Verify the function returns the correct types
        assert isinstance(D, float)
        assert isinstance(p, float)

    # Test with identical distributions (same mean and std) returns low D statistic and high p-value
    def test_identical_distributions_return_low_D_and_high_p(self):
        from src.magnipore import ks_test
        # Set random seed for reproducibility
        np.random.seed(42)

        # Define two identical distributions
        dist1 = (0, 1)  # mean=0, std=1
        dist2 = (0, 1)  # mean=0, std=1

        # Call the function
        D, p = ks_test(dist1, dist2)

        # With identical distributions and seed, we expect a low D value (close to 0)
        # and a high p-value (indicating distributions are similar)
        assert 0 <= D < 0.5, f"Expected low D statistic, got {D}"
        assert p > 0.5, f"Expected high p-value, got {p}"

        # Verify the function returns the correct types
        assert isinstance(D, float)
        assert isinstance(p, float)

    # Calculate td-score for positive mDiff and positive sAvg
    def test_positive_mdiff_and_savg(self):
        from src.magnipore import td_score
        # Arrange
        mdiff = 10.0
        savg = 5.0
        expected_score = 2.0
    
        # Act
        result = td_score(mdiff, savg)
    
        # Assert
        assert result == expected_score

    # Returns correct KL divergence when both standard deviations are positive
    def test_correct_kl_divergence_with_positive_std_devs(self):
        from src.magnipore import kullback_leibler_normal
        # Test with known values
        m0, s0 = 0.0, 1.0
        m1, s1 = 2.0, 2.0
    
        # Calculate expected KL divergence manually
        ratio = (s0 / s1) ** 2
        expected = (ratio + (m1-m0) ** 2/s1 ** 2 - 1 + np.log(ratio)) / 2
    
        # Get actual result
        result = kullback_leibler_normal(m0, s0, m1, s1)
    
        # Assert that the result matches the expected value
        np.testing.assert_almost_equal(result, expected)
    
        # Test with another set of values
        m0, s0 = 5.0, 2.5
        m1, s1 = 7.0, 1.5
    
        # Calculate expected KL divergence manually
        ratio = (s0 / s1) ** 2
        expected = (ratio + (m1-m0) ** 2/s1 ** 2 - 1 + np.log(ratio)) / 2
    
        # Get actual result
        result = kullback_leibler_normal(m0, s0, m1, s1)
    
        # Assert that the result matches the expected value
        np.testing.assert_almost_equal(result, expected)

    # Correctly calculates and compares signal metrics (td score, KL divergence, Bayesian p) for two data positions using a lock that supports context manager protocol.
    def test_signal_metrics_calculation_with_context_manager_lock(self, mocker):
        from src.magnipore import kullback_leibler_normal
        from statistics import NormalDist
        # Mock data positions
        data_pos1 = mocker.Mock()
        data_pos1.get_mean_stdev.return_value = (10.0, 2.0)
        data_pos1.n_reads = 20
        data_pos1.magnipore_string.return_value = "pos1_magnipore_string"

        data_pos2 = mocker.Mock()
        data_pos2.get_mean_stdev.return_value = (12.0, 3.0)
        data_pos2.n_reads = 15
        data_pos2.magnipore_string.return_value = "pos2_magnipore_string"

        # Mock multiprocessing objects
        lock = mocker.MagicMock()  # Use MagicMock to support context manager
        all_queue = mocker.Mock()
        sign_queue = mocker.Mock()
        stk_queue = mocker.Mock()

        # Mock shared values
        no_data = mocker.Mock()
        no_data.value = 0
        low_cov_count = mocker.Mock()
        low_cov_count.value = 0
        num_muts = mocker.Mock()
        num_muts.value = 0
        sign_pos = mocker.Mock()
        sign_pos.value = 0
        num_pos = mocker.Mock()
        num_pos.value = 0

        # Prepare arguments
        args = (
            0,  # strand (+ strand)
            "A", "G",  # base1, base2
            "AATCG", "GATCG",  # motif1, motif2
            "alignment_info",  # alip
            data_pos1, data_pos2,  # data positions
            ["sample1", "sample2"],  # seqs_ids
            100, 200,  # pos1, pos2
            num_muts, sign_pos, no_data, low_cov_count, num_pos,
            lock, all_queue, sign_queue, stk_queue
        )

        # Call the function
        from src.magnipore import compare_signals
        compare_signals(args)

        # Expected values
        expected_td = abs(10.0 - 12.0) / ((2.0 + 3.0) / 2)  # 0.8
        expected_kl = kullback_leibler_normal(10.0, 2.0, 12.0, 3.0)
        expected_bayesian_p = NormalDist(10.0, 2.0).overlap(NormalDist(12.0, 3.0))

        # Verify queue put was called with correct data
        all_queue.put.assert_called_once()
        put_args = all_queue.put.call_args[0][0]

        # Check that metrics are in the output string
        assert f"{expected_td:.8f}" in put_args
        assert f"{expected_kl:.8f}" in put_args
        assert f"{expected_bayesian_p:.8f}" in put_args

        # Verify sign_queue was not called (td < 1)
        sign_queue.put.assert_not_called()

        # Verify counters were updated correctly
        assert num_pos.value == 1
        assert low_cov_count.value == 0  # Both have coverage > 10

    # Parse command line arguments with all required arguments provided
    def test_parse_with_all_required_arguments(self):
        from unittest.mock import patch
        from src.magnipore import parse
    
        # Setup command line arguments
        test_args = [
            'magnipore',
            '/path/to/raw_data_first.pod5',
            '/path/to/raw_data_second.pod5',
            '/path/to/basecalls_first.bam',
            '/path/to/basecalls_second.bam',
            '/path/to/uncalled4_first',
            '/path/to/uncalled4_second',
            '/path/to/alignment.bam',
            '/path/to/outdir',
            'dna_r9'
        ]
    
        with patch('sys.argv', test_args):
            args = parse()
        
            # Assert all required arguments are correctly parsed
            assert args.raw_data_first_sample == '/path/to/raw_data_first.pod5'
            assert args.raw_data_sec_sample == '/path/to/raw_data_second.pod5'
            assert args.basecalls_data_first_sample == '/path/to/basecalls_first.bam'
            assert args.basecalls_data_sec_sample == '/path/to/basecalls_second.bam'
            assert args.uncalled4_first_sample == '/path/to/uncalled4_first'
            assert args.uncalled4_sec_sample == '/path/to/uncalled4_second'
            assert args.alignment == '/path/to/alignment.bam'
            assert args.outdir == '/path/to/outdir'
            assert args.pore == 'dna_r9'
        
            # Assert default values for optional arguments
            assert args.label_first_sample == 'sample_1'
            assert args.label_sec_sample == 'sample_2'
            assert args.threads == 1
            assert args.calculate_data_density is False

class TestHelper:
    # Format a small number of bytes (e.g., 100.0) with default suffix
    def test_small_number_with_default_suffix(self):
        from src.Helper import sizeof_fmt
        # Test that a small number of bytes is formatted correctly with default suffix
        result = sizeof_fmt(100.0)
        assert result == "100.0B"

        # Format a very large number that reaches 'Yi' prefix (>= 1024^8)
    def test_very_large_number_reaches_yi_prefix(self):
        from src.Helper import sizeof_fmt
        # Test that a very large number (>= 1024^8) is formatted with 'Yi' prefix
        # 1024^8 = 1,180,591,620,717,411,303,424
        large_number = 1024**8
        result = sizeof_fmt(large_number)
        assert result == "1.0YiB"

    # Returns the correct reverse complement of a DNA sequence based on complementary base pairing rules.
    def test_correct_reverse_complement(self):
        from src.Helper import rev_complement

        # Test with standard DNA sequence
        input_seq = "ACGT"
        expected = "ACGT"[::-1].translate(str.maketrans("ACGT", "TGCA"))
        result = rev_complement(input_seq)
        assert result == expected

        # Test with longer standard DNA sequence
        input_seq = "AAAAACCCCCGGGGGTTTTTT"
        expected = "AAAAACCCCCGGGGGTTTTTT"[::-1].translate(str.maketrans("ACGT", "TGCA"))
        result = rev_complement(input_seq)
        assert result == expected

    # Correctly handles non-standard nucleotide codes (Y, R, S, W, M, K) and returns the expected reverse complement.
    def test_nonstandard_nucleotide_codes_fixed(self):
        from src.Helper import rev_complement

        # Test with non-standard nucleotide codes
        input_seq = "YRSWMK"
        expected = "MKWSYR"
        result = rev_complement(input_seq)

        assert result == expected

        # Test with mixed standard and non-standard codes
        input_seq = "ACGTYRMK"
        expected = "MKYRACGT"
        result = rev_complement(input_seq)

        assert result == expected

    # Complementing a standard DNA sequence with A, C, G, T returns correct complements
    def test_standard_dna_sequence_complementation(self):
        from src.Helper import complement

        # Test with standard DNA sequence
        input_seq = "ACGT"
        expected_output = "TGCA"

        result = complement(input_seq)

        assert result == expected_output

    # Complementing a sequence with characters not in COMPLEMENT dictionary returns 'N' for those characters
    def test_unknown_characters_return_n(self):
        from src.Helper import complement

        # Test with sequence containing characters not in COMPLEMENT dictionary
        input_seq = "ACGTXZ123"
        expected_output = "TGCANNNNN"

        result = complement(input_seq)

        assert result == expected_output

class TestRed:
    # Initializing a Red object with default parameters
    def test_default_initialization(self):
        import numpy as np
        from src.Red import Red

        # Initialize Red with default parameters
        red = Red()

        # Check default values
        assert red.INITLEN == 30
        assert red.SKIP_CALC is False
        assert red.k == 0
        assert red.n == 0
        assert red.ex == 0.0
        assert red.ex2 == 0.0
        assert np.array_equal(red.initvec, np.zeros(0, dtype=np.float32))
        assert red.data_density == 0.0
        assert red.n_datapoints == 0
        assert red.contained_datapoints == 0
        assert red.n_segments == 0
        assert red.contained_segments == 0
        assert red.n_reads == 0
        assert red.mean is None
        assert red.std is None
        assert red.var is None

    # Verify that appending an empty array results in a NaN mean and zero standard deviation
    @pytest.mark.filterwarnings("ignore:Mean of empty slice")
    def test_append_empty_array_with_nan_mean(self):
        import numpy as np
        from src.Red import Red

        # Initialize Red
        red = Red()

        # Append an empty array
        empty_array = np.array([], dtype=np.float32)
        red.append(empty_array)

        # Check that nothing changed in the statistics
        assert red.n == 0
        assert red.ex == 0.0
        assert red.ex2 == 0.0

        # Check that the empty array was appended to initvec
        assert np.array_equal(red.initvec, empty_array)

        # Verify that _flush wasn't triggered (since buffer isn't full)
        assert len(red.initvec) == 0

        # Get mean and stdev should handle this gracefully
        mean, std = red.get_mean_stdev()
        assert np.isnan(mean)  # Mean should be NaN when no data is processed
        assert std == 0.0  # Default when n < 2

    # Appending values to the model and calculating mean/standard deviation
    def test_append_and_calculate_mean_stdev(self):
        import numpy as np
        from src.Red import Red

        # Initialize Red with default parameters
        red = Red()

        # Append values to the model
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
        red.append(values)

        # Calculate mean and standard deviation
        mean, std = red.get_mean_stdev()

        # Check if the mean and standard deviation are calculated correctly
        expected_mean = np.mean(values)
        expected_std = np.std(values, ddof=1)  # Sample standard deviation

        assert np.isclose(mean, expected_mean), f"Expected mean: {expected_mean}, but got: {mean}"
        assert np.isclose(std, expected_std), f"Expected std: {expected_std}, but got: {std}"

    # Tracking read, segment, and datapoint counts through add methods
    def test_add_methods_update_counts(self):
        from src.Red import Red
        # Initialize Red with default parameters
        red = Red()

        # Add reads, segments, and datapoints
        red.add_reads(5)
        red.add_segments(3)
        red.add_datapoints(10)
        red.add_contained_datapoints(7)
        red.add_contained_segments(2)

        # Check if the counts are updated correctly
        assert red.n_reads == 5
        assert red.n_segments == 3
        assert red.n_datapoints == 10
        assert red.contained_datapoints == 7
        assert red.contained_segments == 2

    # Setting and calculating data density values
    def test_set_and_add_data_density(self):
        from src.Red import Red
        

        # Initialize Red with default parameters
        red = Red()

        # Set data density and verify
        red.set_data_density(5.0)
        assert red.data_density == 5.0

        # Add to data density and verify
        red.add_data_density(2.5)
        assert red.data_density == 7.5

    # Getting string representations for debugging and Magnipore output
    def test_string_representations(self):
        import numpy as np
        from src.Red import Red

        # Initialize Red with default parameters
        red = Red()

        # Append some data to trigger calculations
        red.append(np.array([1.0, 2.0, 3.0]))
        red.append(np.array([4.0, 5.0, 6.0]))

        # Get string representations
        debug_str = str(red)
        magnipore_str = red.magnipore_string()

        # Check if the string representations are as expected
        assert isinstance(debug_str, str)
        assert isinstance(magnipore_str, str)

        # Check if the debug string contains expected values
        assert f"{red.mean:.8f}" in debug_str
        assert f"{red.std:.8f}" in debug_str
        assert f"{red.data_density:.8f}" in debug_str

        # Check if the magnipore string contains expected values
        assert f"{red.n_datapoints:.0f}" in magnipore_str
        assert f"{red.contained_datapoints:.0f}" in magnipore_str
        assert f"{red.n_segments:.0f}" in magnipore_str
        assert f"{red.contained_segments:.0f}" in magnipore_str
        assert f"{red.n_reads:.0f}" in magnipore_str

    # Online mean-variance tracking with multiple data batches
    def test_online_mean_variance_tracking_multiple_batches(self):
        import numpy as np
        from src.Red import Red

        # Initialize Red with default parameters
        red = Red()

        # Append a batch of numbers
        data1 = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0,
                            11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
                            21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0])
        red.append(data1)

        # Force flush to process the data
        red._flush()

        # Calculate expected mean and variance for the first batch
        expected_mean1 = np.mean(data1)
        expected_variance1 = np.var(data1, ddof=1)  # Sample variance

        # Retrieve calculated mean and standard deviation after first batch
        mean1, std1 = red.get_mean_stdev()

        # Assert the mean and variance are as expected for the first batch
        assert np.isclose(mean1, expected_mean1), f"Expected mean: {expected_mean1}, but got: {mean1}"
        assert np.isclose(std1, np.sqrt(expected_variance1)), f"Expected std: {np.sqrt(expected_variance1)}, but got: {std1}"

        # Append another batch of numbers
        data2 = np.array([31.0, 32.0, 33.0, 34.0, 35.0])
        red.append(data2)

        # Force flush to process the second batch of data
        red._flush()

        # Combine data for expected calculations
        combined_data = np.concatenate((data1, data2))

        # Calculate expected mean and variance for combined data
        expected_mean_combined = np.mean(combined_data)
        expected_variance_combined = np.var(combined_data, ddof=1)  # Sample variance

        # Retrieve calculated mean and standard deviation after second batch
        mean_combined, std_combined = red.get_mean_stdev()

        # Assert the mean and variance are as expected for combined data
        assert np.isclose(mean_combined, expected_mean_combined), f"Expected mean: {expected_mean_combined}, but got: {mean_combined}"
        assert np.isclose(std_combined, np.sqrt(expected_variance_combined)), f"Expected std: {np.sqrt(expected_variance_combined)}, but got: {std_combined}"

    # Online mean-variance tracking with multiple batches of random numbers
    def test_online_mean_variance_tracking_random_batches(self):
        import numpy as np
        from src.Red import Red

        # Initialize Red with default parameters
        red = Red()

        # Append a batch of random numbers
        data1 = np.random.rand(30) * 100  # Random numbers between 0 and 100
        red.append(data1)

        # Force flush to process the data
        red._flush()

        # Calculate expected mean and variance for the first batch
        expected_mean1 = np.mean(data1)
        expected_variance1 = np.var(data1, ddof=1)  # Sample variance

        # Retrieve calculated mean and standard deviation after first batch
        mean1, std1 = red.get_mean_stdev()

        # Assert the mean and variance are as expected for the first batch
        assert np.isclose(mean1, expected_mean1), f"Expected mean: {expected_mean1}, but got: {mean1}"
        assert np.isclose(std1, np.sqrt(expected_variance1)), f"Expected std: {np.sqrt(expected_variance1)}, but got: {std1}"

        # Append another batch of random numbers
        data2 = np.random.rand(5) * 100  # Random numbers between 0 and 100
        red.append(data2)

        # Force flush to process the second batch of data
        red._flush()

        # Combine data for expected calculations
        combined_data = np.concatenate((data1, data2))

        # Calculate expected mean and variance for combined data
        expected_mean_combined = np.mean(combined_data)
        expected_variance_combined = np.var(combined_data, ddof=1)  # Sample variance

        # Retrieve calculated mean and standard deviation after second batch
        mean_combined, std_combined = red.get_mean_stdev()

        # Assert the mean and variance are as expected for combined data
        assert np.isclose(mean_combined, expected_mean_combined), f"Expected mean: {expected_mean_combined}, but got: {mean_combined}"
        assert np.isclose(std_combined, np.sqrt(expected_variance_combined)), f"Expected std: {np.sqrt(expected_variance_combined)}, but got: {std_combined}"

    # Testing buffer management with values less than initlen
    def test_append_with_values_less_than_initlen(self):
        import numpy as np
        from src.Red import Red
    
        # Initialize Red with default parameters
        red = Red(initlen=5)
    
        # Append values less than initlen
        red.append(np.array([1.0, 2.0, 3.0]))
    
        # Check that the buffer contains the appended values
        assert np.array_equal(red.initvec, np.array([1.0, 2.0, 3.0], dtype=np.float32))
    
        # Ensure that _flush() has not been called
        assert red.n == 0
        assert red.ex == 0.0
        assert red.ex2 == 0.0

    # Avoids unnecessary recalculation when SKIP_CALC is True
    def test_avoids_recalculation_when_skip_calc_is_true(self):
        from src.Red import Red
        # Arrange
        red = Red(skip_calc=True)
        red.mean = 3.0
        red.std = 1.5
    
        # Act
        actual_mean, actual_std = red.get_mean_stdev()
    
        # Assert
        assert actual_mean == 3.0
        assert actual_std == 1.5

    # Returns samples from the reservoir when the reservoir contains samples
    def test_returns_samples_when_reservoir_contains_data(self, mocker):
        # Arrange
        mock_reservoir = mocker.Mock()
        mock_samples = np.array([1, 2, 3, 4, 5])
        mock_reservoir.samples.return_value = mock_samples
    
        from src.Red import Red
        red = Red()
        red.reservoir = mock_reservoir
    
        # Act
        result = red.get_samples()
    
        # Assert
        assert np.array_equal(result, mock_samples)
        mock_reservoir.samples.assert_called_once()

    # Returns empty array when reservoir is empty
    def test_returns_empty_array_when_reservoir_is_empty(self, mocker):
        # Arrange
        mock_reservoir = mocker.Mock()
        mock_reservoir.samples.return_value = np.array([])
    
        from src.Red import Red
        red = Red()
        red.reservoir = mock_reservoir
    
        # Act
        result = red.get_samples()
    
        # Assert
        assert isinstance(result, np.ndarray)
        assert result.size == 0
        mock_reservoir.samples.assert_called_once()

class TestLogger:
    # Logger initializes with a logfilepointer and s logs to it
    def test_logger_with_logfilepointer_writes_logs(self, mocker):
        from src.Logger import Logger
        # Arrange
        mock_file = mocker.mock_open()
        mock_file_instance = mock_file.return_value
        logger = Logger(mock_file_instance)

        # Act
        logger._writeLog("Test log message")

        # Assert
        mock_file_instance.write.assert_called_once_with("Test log message")

    # Logger.error() writes to stderr and logfile with timestamp and mocked memory usage
    def test_error_writes_to_stderr_and_logfile(self, mocker):
        from src.Logger import Logger
        # Arrange
        mock_file = mocker.mock_open()
        mock_file_instance = mock_file.return_value
        logger = Logger(mock_file_instance)
        mock_stderr = mocker.patch('sys.stderr', new_callable=mocker.Mock())
        mock_exit = mocker.patch('sys.exit')
        mock_datetime = mocker.patch('datetime.datetime')
        mock_datetime.now.return_value.strftime.return_value = '2023-10-01_12-00-00'
        # Correct the mock path for memory usage
        mocker.patch('psutil.Process.memory_info', return_value=mocker.Mock(rss=100 * 1024 * 1024))

        # Act
        logger.error("Test error message", error_type='1')

        # Assert
        expected_output = '2023-10-01_12-00-00, MEM: 100.00 MB, \033[91mERROR: Test error message\nMagnipore Error Code: 1\n\033[0m\n'
        mock_stderr.write.assert_called_once_with(expected_output)
        mock_file_instance.write.assert_called_once_with('2023-10-01_12-00-00, MEM: 100.00 MB, ERROR: Test error message\nMagnipore Error Code: 1\n')
        mock_exit.assert_called_once_with('1')

    # Logger.warning() writes to stderr and logfile with timestamp and memory usage
    def test_warning_writes_to_stderr_and_logfile(self, mocker):
        from src.Logger import Logger, get_memory_usage
        import datetime
        from src.Helper import ANSI
        
        # Arrange
        mock_file = mocker.mock_open()
        mock_file_instance = mock_file.return_value
        logger = Logger(mock_file_instance)
        mock_stderr = mocker.patch('sys.stderr', new_callable=mocker.Mock())
        test_message = "Test warning message"
        timestamp = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
        memory_usage = f"{get_memory_usage():.2f} MB"
    
        # Act
        logger.warning(test_message)
    
        # Assert
        expected_output = f'{timestamp}, MEM: {memory_usage}, {ANSI.RED}WARNING: {test_message}{ANSI.END}\n'
        mock_stderr.write.assert_called_once_with(expected_output)
        mock_file_instance.write.assert_called_once_with(f'{timestamp}, MEM: {memory_usage}, WARNING: {test_message}\n')

    # Logger.writeLog() writes string to logfile when logfilepointer is set
    def test_write_log_writes_to_logfile_when_pointer_is_set(self, mocker):
        from src.Logger import Logger
        # Arrange
        mock_file = mocker.mock_open()
        mock_file_instance = mock_file.return_value
        logger = Logger(mock_file_instance)
    
        # Act
        logger._writeLog("Test log message")
    
        # Assert
        mock_file_instance.write.assert_called_once_with("Test log message")

    # Logger.printLog() writes to stdout and logfile with timestamp and memory usage
    def test_printlog_writes_to_stdout_and_logfile(self, mocker):
        from src.Logger import Logger, get_memory_usage
        import datetime
        # Arrange
        mock_file = mocker.mock_open()
        mock_file_instance = mock_file.return_value
        logger = Logger(mock_file_instance)
        mock_stdout_write = mocker.patch('sys.stdout.write')
        test_string = "Test log message"

        # Act
        logger.printLog(test_string)

        # Assert
        expected_output = f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}, MEM: {get_memory_usage():.2f} MB, LOG: {test_string}\n'
        mock_stdout_write.assert_any_call(expected_output.strip())
        mock_file_instance.write.assert_any_call(expected_output.strip())

    # Logger.printLog() correctly writes to stdout and logfile with newline_before=True and newline_after=False
    def test_print_log_with_newline_before_and_no_newline_after(self, mocker):
        from src.Logger import Logger, get_memory_usage
        import datetime
        # Arrange
        mock_stdout_write = mocker.patch('sys.stdout.write')
        mock_file = mocker.mock_open()
        mock_file_instance = mock_file.return_value
        logger = Logger(mock_file_instance)
        test_string = "Test log message"

        # Act
        logger.printLog(test_string, newline_before=True, newline_after=False)

        # Assert
        expected_log_output = f'{datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")}, MEM: {get_memory_usage():.2f} MB, LOG: {test_string}'
        mock_stdout_write.assert_any_call('\n')
        mock_stdout_write.assert_any_call(expected_log_output)
        mock_file_instance.write.assert_any_call('\n')
        mock_file_instance.write.assert_any_call(expected_log_output)

    # Returns memory usage in MB for the current process
    def test_returns_memory_usage_in_mb(self, mocker):
        from src.Logger import get_memory_usage
        # Arrange
        mock_process = mocker.Mock()
        mock_process.memory_info.return_value.rss = 104857600  # 100 MB in bytes
        mocker.patch('psutil.Process', return_value=mock_process)
    
        # Act
        result = get_memory_usage()
    
        # Assert
        assert result == 100.0  # Should return 100 MB
        mock_process.memory_info.assert_called_once()

    # Handles large memory usage values (several GB)
    def test_handles_large_memory_values(self, mocker):
        from src.Logger import get_memory_usage
        # Arrange
        mock_process = mocker.Mock()
        mock_process.memory_info.return_value.rss = 5 * 1024 * 1024 * 1024  # 5 GB in bytes
        mocker.patch('psutil.Process', return_value=mock_process)
    
        # Act
        result = get_memory_usage()
    
        # Assert
        assert result == 5120.0  # Should return 5120 MB (5 GB)
        mock_process.memory_info.assert_called_once()

class TestReservoir:
    
    # Adding elements to an empty reservoir until it's full
    def test_adding_elements_until_full(self):
        import numpy as np
        from src.Reservoir import Reservoir
    
        # Initialize reservoir with k=5
        k = 5
        reservoir = Reservoir(k)
    
        # Add elements one by one until full
        for i in range(k):
            reservoir.add(np.array([i]))
        
            # Check that the reservoir contains the correct elements
            samples = reservoir.samples()
            assert len(samples) == i + 1
            assert np.array_equal(samples, np.array(range(i + 1)))
        
        # Verify the reservoir is now full
        assert reservoir.cnt == k
        assert len(reservoir.samples()) == k
        assert np.array_equal(reservoir.samples(), np.array(range(k)))

    # Adding elements to a reservoir that's already full
    def test_adding_elements_to_full_reservoir(self):
        import numpy as np
        from src.Reservoir import Reservoir
    
        # Initialize reservoir with k=5
        k = 5
        reservoir = Reservoir(k)
    
        # Fill the reservoir
        for i in range(k):
            reservoir.add(np.array([i]))
    
        # Add additional elements beyond the capacity
        additional_elements = np.array([5, 6, 7, 8, 9])
        reservoir.add(additional_elements)
    
        # Check that the reservoir still contains k elements
        samples = reservoir.samples()
        assert len(samples) == k
    
        # Since the reservoir is full, we expect the elements to be probabilistically replaced
        # We cannot predict the exact content, but we can check that all elements are within the range of added elements
        for sample in samples:
            assert sample in np.concatenate((np.array(range(k)), additional_elements))

    # Retrieving samples from a partially filled reservoir
    def test_retrieving_samples_from_partially_filled_reservoir(self):
        import numpy as np
        from src.Reservoir import Reservoir
    
        # Initialize reservoir with k=5
        k = 5
        reservoir = Reservoir(k)
    
        # Add fewer elements than k
        elements_to_add = 3
        for i in range(elements_to_add):
            reservoir.add(np.array([i]))
    
        # Retrieve samples from the partially filled reservoir
        samples = reservoir.samples()
    
        # Check that the reservoir contains the correct number of elements
        assert len(samples) == elements_to_add
        assert np.array_equal(samples, np.array(range(elements_to_add)))

    # Adding empty arrays to the reservoir
    def test_adding_empty_arrays(self):
        import numpy as np
        from src.Reservoir import Reservoir
    
        # Initialize reservoir with k=5
        k = 5
        reservoir = Reservoir(k)
    
        # Add an empty array
        reservoir.add(np.array([]))
    
        # Check that the reservoir is still empty
        samples = reservoir.samples()
        assert len(samples) == 0
        assert np.array_equal(samples, np.array([]))
    
        # Verify the internal counter is still zero
        assert reservoir.cnt == 0

    # Checking if the reservoir maintains exactly k samples after many additions
    def test_reservoir_maintains_k_samples_after_many_additions(self):
        import numpy as np
        from src.Reservoir import Reservoir
    
        # Initialize reservoir with k=5
        k = 5
        reservoir = Reservoir(k)
    
        # Add more than k elements
        num_elements = 100
        elements = np.arange(num_elements)
        reservoir.add(elements)
    
        # Check that the reservoir contains exactly k samples
        samples = reservoir.samples()
        assert len(samples) == k
        assert all(sample in elements for sample in samples)

    # Retrieving samples from a completely filled reservoir
    def test_retrieving_samples_from_filled_reservoir(self):
        import numpy as np
        from src.Reservoir import Reservoir
    
        # Initialize reservoir with k=5
        k = 5
        reservoir = Reservoir(k)
    
        # Add more than k elements to fill the reservoir
        elements = np.array(range(10))
        reservoir.add(elements)
    
        # Retrieve samples from the filled reservoir
        samples = reservoir.samples()
    
        # Verify the reservoir is full and contains k samples
        assert len(samples) == k
        assert all(sample in elements for sample in samples)

    # Verifying the randomness of the sampling algorithm
    def test_randomness_of_sampling(self):
        import numpy as np
        from src.Reservoir import Reservoir
    
        # Initialize reservoir with k=5
        k = 5
        reservoir = Reservoir(k)
    
        # Add more elements than the reservoir can hold
        num_elements = 1000
        elements = np.arange(num_elements)
        reservoir.add(elements)
    
        # Retrieve samples from the reservoir
        samples = reservoir.samples()
    
        # Check that the reservoir contains exactly k elements
        assert len(samples) == k
    
        # Check that the samples are a subset of the added elements
        assert all(sample in elements for sample in samples)
    
        # Check that the samples are not in a predictable order
        # This is a simple check to ensure randomness, not a statistical test
        assert not np.array_equal(samples, np.arange(k))

    # Testing if the reservoir sampling is unbiased using a chi-square test
    # def test_reservoir_sampling_unbiasedness_with_chi_square(self):
    #     import numpy as np
    #     from collections import Counter
    #     from scipy.stats import chisquare
    #     from src.Reservoir import Reservoir

    #     # Initialize parameters
    #     k = 10
    #     n = 1000
    #     trials = 10000

    #     # Create a large stream of elements
    #     stream = np.arange(n)

    #     # Counter to track occurrences of each element in the reservoir
    #     element_counts = Counter()

    #     # Perform multiple trials to check unbiasedness
    #     for _ in range(trials):
    #         reservoir = Reservoir(k)
    #         reservoir.add(stream)
    #         samples = reservoir.samples()
    #         element_counts.update(samples)

    #     # Calculate expected count for each element
    #     expected_count = [trials * k / n] * n

    #     # Perform chi-square test
    #     observed_counts = [element_counts[i] for i in range(n)]
    #     chi2, p_value = chisquare(observed_counts, expected_count)

    #     # Assert that the p-value is greater than 0.05 for unbiasedness
    #     assert p_value > 0.05

    # Verifying the internal state counters (cnt, next) are correctly updated
    def test_internal_counters_update_correctly(self):
        import numpy as np
        from src.Reservoir import Reservoir
    
        # Initialize reservoir with k=5
        k = 5
        reservoir = Reservoir(k)
    
        # Add elements one by one and check internal counters
        for i in range(10):
            reservoir.add(np.array([i]))
        
            # Check that the counter is incremented correctly
            assert reservoir.cnt == i + 1
        
            # Check that 'next' is updated correctly after filling the reservoir
            if reservoir.cnt >= k:
                assert reservoir.next > reservoir.cnt

class TestMagniplot:

    # Ensure the marker function returns 'D' for uppercase and mixed case mutation contexts.
    def test_marker(self):
        from src.magniplot import marker
        # Act
        result_upper = marker('MUT')
        result_mixed = marker('Mut')

        # Assert
        assert result_upper == 'D', "Should return 'D' for 'MUT' (case insensitive)"
        assert result_mixed == 'D', "Should return 'D' for 'Mut' (case insensitive)"

        # Act
        result_upper = marker('MOD')
        result_mixed = marker('mod')

        # Assert
        assert result_upper == 'o', "Should return 'o' for 'MOD' (case insensitive)"
        assert result_mixed == 'o', "Should return 'o' for 'mod' (case insensitive)"

    # Returns 'blue' when mut_context is 'mut'
    def test_color(self):
        from src.magniplot import color
        # Act
        result_upper = color('MUT')
        result_mixed = color('Mut')

        # Assert
        assert result_upper == 'blue', "Should return 'blue' for 'MUT' (case insensitive)"
        assert result_mixed == 'blue', "Should return 'blue' for 'Mut' (case insensitive)"

        # Act
        result_upper = color('MOD')
        result_mixed = color('mod')

        # Assert
        assert result_upper == '#d95f02', "Should return '#d95f02' for 'MOD' (case insensitive)"
        assert result_mixed == '#d95f02', "Should return '#d95f02' for 'mod' (case insensitive)"