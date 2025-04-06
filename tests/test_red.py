import pytest
from src.Red import Red
from src.OnlineMeanVar import OnlineMeanVar
import numpy as np

class TestRed:

    # Initializing Red with default parameters creates an instance with correct initial values
    def test_init_default_parameters(self):
        # Arrange & Act
        red = Red()
    
        # Assert
        assert red.data_density == 0.0
        assert red.n_datapoints == 0
        assert red.contained_datapoints == 0
        assert red.n_segments == 0
        assert red.contained_segments == 0
        assert red.n_reads == 0
    
        # Check OnlineMeanVar instances were created with correct parameters
        assert isinstance(red.signal_stats, OnlineMeanVar)
        assert red.signal_stats.INITLEN == 30
        assert red.signal_stats.SKIP_CALC is False
    
        assert isinstance(red.dwell_time_stats, OnlineMeanVar)
        assert red.dwell_time_stats.INITLEN == 30
        assert red.dwell_time_stats.SKIP_CALC is False

    # Initializing with skip_calc=True should bypass calculations in get_mean_stdev
    def test_skip_calc_bypasses_calculations(self):
        # Arrange
        red = Red(skip_calc=True)
        test_data = np.array([1.0, 2.0, 3.0], dtype=np.float32)
    
        # Act
        red.append(test_data)
    
        # Get the original mean and std values
        red.signal_stats.mean = 10.0
        red.signal_stats.std = 5.0
    
        # Call get_mean_stdev which should return the preset values without recalculation
        mean, std = red.get_signal_mean_stdev()
    
        # Assert
        assert mean == 10.0
        assert std == 5.0
    
        # Verify the same behavior for dwell_time_stats
        red.dwell_time_stats.mean = 20.0
        red.dwell_time_stats.std = 7.0
    
        dwell_mean, dwell_std = red.get_dwell_time_mean_stdev()
        assert dwell_mean == 20.0
        assert dwell_std == 7.0

    # Adding reads, segments, datapoints increments respective counters by specified amounts
    def test_increment_counters(self):
        # Arrange
        red = Red()
    
        # Act
        red.add_reads(5)
        red.add_segments(3)
        red.add_datapoints(10)
        red.add_contained_datapoints(7)
        red.add_contained_segments(2)
    
        # Assert
        assert red.n_reads == 5
        assert red.n_segments == 3
        assert red.n_datapoints == 10
        assert red.contained_datapoints == 7
        assert red.contained_segments == 2

    # magnipore_string returns formatted string with count statistics
    def test_magnipore_string_format(self):
        # Arrange
        red = Red()
        red.add_datapoints(100)
        red.add_contained_datapoints(80)
        red.add_segments(10)
        red.add_contained_segments(8)
        red.add_reads(5)
    
        # Act
        result = red.magnipore_string()
    
        # Assert
        expected_output = "100\t80\t10\t8\t5"
        assert result == expected_output

    # Returns a formatted string with all statistics when all values are valid
    def test_returns_formatted_string_with_valid_values(self, mocker):
        # Arrange
        from src.Red import Red
        from src.OnlineMeanVar import OnlineMeanVar
    
        # Mock OnlineMeanVar instances
        mock_signal_stats = mocker.Mock(spec=OnlineMeanVar)
        mock_signal_stats.get_mean_stdev.return_value = (1.234, 0.567)
    
        mock_dwell_time_stats = mocker.Mock(spec=OnlineMeanVar)
        mock_dwell_time_stats.get_mean_stdev.return_value = (10.987, 2.345)
    
        # Create Red instance with mocked dependencies
        red = Red()
        red.signal_stats = mock_signal_stats
        red.dwell_time_stats = mock_dwell_time_stats
    
        # Set other required attributes
        red.data_density = 0.789
        mocker.patch.object(red, 'expected_model_density', return_value=0.321)
        red.n_datapoints = 100
        red.contained_datapoints = 90
        red.n_segments = 5
        red.contained_segments = 4
        red.n_reads = 10
    
        # Act
        result = str(red)
    
        # Assert
        expected = (
            f"1.23400000\t0.56700000\t"
            f"10.98700000\t2.34500000\t"
            f"0.78900000\t0.32100000\t"
            f"100\t90\t"
            f"5\t4\t"
            f"10"
        )
        assert result == expected