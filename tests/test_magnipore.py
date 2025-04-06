import pytest
import multiprocessing
import numpy as np

multiprocessing.set_start_method("forkserver", force=True)  # Alternative: "spawn"

class TestMagnipore:
    @pytest.fixture(scope="module")
    def test_data(self):
        from Bio import SeqIO
        import pandas as pd
        from src.magnipore import init_Logger, read_red_file, getMapping
        from os.path import join, dirname

        """Load test data from files and return a dictionary of necessary objects."""
        init_Logger(None)
        this_file_dir = dirname(__file__)

        # File paths
        ref1 = join(this_file_dir, 'sample1.fa')
        red1 = join(this_file_dir, 'sample1.red')
        ref2 = join(this_file_dir, 'sample2.fa')
        red2 = join(this_file_dir, 'sample2.red')
        alignment_path = join(this_file_dir, 'alignment', 'sample1_sample2.aln')

        # Read sequences and red files
        seq_dict = {**SeqIO.to_dict(SeqIO.parse(ref1, "fasta")), **SeqIO.to_dict(SeqIO.parse(ref2, "fasta"))}
        red1DF = pd.read_csv(red1, sep='\t')
        red2DF = pd.read_csv(red2, sep='\t')

        # Process alignment
        mapping_dict, unaligned_dict, aln_dict, seq_dict = getMapping(alignment_path, this_file_dir, 'sample1', 'sample2')
        red1_dict = read_red_file(red1, list(seq_dict.values())[0])
        red2_dict = read_red_file(red2, list(seq_dict.values())[-1])

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
                assert np.isclose(red.signal_stats.mean, red1DF.loc[(2*i)+j, 'signal_mean']), "Mismatch in signal_mean"
                assert np.isclose(red.signal_stats.std, red1DF.loc[(2*i)+j, 'signal_std']), "Mismatch in signal_std"
                assert np.isclose(red.dwell_time_stats.mean, red1DF.loc[(2*i)+j, 'dwell_time_mean']), "Mismatch in signal_mean"
                assert np.isclose(red.dwell_time_stats.std, red1DF.loc[(2*i)+j, 'dwell_time_std']), "Mismatch in signal_std"
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
        from src.magnipore import magnipore
        from os.path import join
        from pandas import read_csv
        """Test the `magnipore` function by validating output file contents."""
        magnipore(
            test_data["mapping_dict"], test_data["unaligned_dict"],
            test_data["seq_dict"], test_data["aln_dict"],
            test_data["red1_dict"], test_data["red2_dict"],
            'sample1', 'sample2', test_data["this_file_dir"], 5, 8
        )

        magn_file = join(test_data["this_file_dir"], 'magnipore/sample1_sample2/sample1_sample2.magnipore')
        magn = read_csv(magn_file, sep='\t')

        assert not magn.empty, "Magnipore output file is empty"
        assert "td_score" in magn.columns, "Missing td_score column in magnipore output"
        assert "bayesian_p" in magn.columns, "Missing bayesian_p column in magnipore output"

    def test_stockholm(self, test_data):
        from os.path import join
        """Test if the marked stockholm file is correctly written with expected base modifications."""
        stk_file = join(test_data["this_file_dir"], 'magnipore/sample1_sample2/sample1_sample2_marked.stk')

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
        
        from src.Red import Red

        # Create original Reds with known values
        original_reds = []
        for i in range(8):
            red = Red()
            # Add different values to each Red
            red.append(np.array([i*10 + j for j in range(5)], dtype=np.float32))
            red.get_signal_mean_stdev()  # Calculate stats
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
            orig_mean, orig_std = red.get_signal_mean_stdev()
            recon_mean, recon_std = reconstructed[i].get_signal_mean_stdev()
    
            assert np.isclose(orig_mean, recon_mean)
            assert np.isclose(orig_std, recon_std)
            assert red.signal_stats.n == reconstructed[i].signal_stats.n, f"Sample count mismatch at index {i}"

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
        red1.__str__.return_value = "10.12345678\t2.12345678\t10.12345678\t2.12345678\t0.12345678\t0.23456789\t100\t80\t5\t4\t10"

        red2 = mocker.MagicMock(spec=Red)
        red2.__str__.return_value = "15.12345678\t3.12345678\t15.12345678\t3.12345678\t0.22345678\t0.33456789\t200\t150\t10\t8\t20"

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
        assert content[0] == 'strand\tposition\tsignal_mean\tsignal_std\tdwell_time_mean\tdwell_time_std\tdata_density\texpected_model_density\tn_datapoints\tcontained_datapoints\tn_segments\tcontained_segments\tn_reads\n'

        # Check data rows
        assert content[1] == '+\t0\t10.12345678\t2.12345678\t10.12345678\t2.12345678\t0.12345678\t0.23456789\t100\t80\t5\t4\t10\n'
        assert content[2] == '+\t1\t15.12345678\t3.12345678\t15.12345678\t3.12345678\t0.22345678\t0.33456789\t200\t150\t10\t8\t20\n'

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
        data_pos1.get_signal_mean_stdev.return_value = (10.0, 2.0)
        data_pos1.n_reads = 20
        data_pos1.magnipore_string.return_value = "pos1_magnipore_string"

        data_pos2 = mocker.Mock()
        data_pos2.get_signal_mean_stdev.return_value = (12.0, 3.0)
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