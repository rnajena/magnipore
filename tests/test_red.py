import pytest
import numpy as np

class TestRed:
    # Initializing a Red object with default parameters
    def test_default_initialization(self):
        
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
