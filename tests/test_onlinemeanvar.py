class TestOnlineMeanVar:

    # Initializing with default parameters creates an empty buffer and zeroed statistics
    def test_default_initialization(self):
        import numpy as np
        from src.OnlineMeanVar import OnlineMeanVar
    
        # Create instance with default parameters
        omv = OnlineMeanVar()
    
        # Check default values
        assert omv.k == 0
        assert omv.n == 0
        assert omv.ex == 0.0
        assert omv.ex2 == 0.0
        assert omv.INITLEN == 30
        assert omv.SKIP_CALC is False
        assert omv.mean is None
        assert omv.std is None
        assert isinstance(omv.initvec, np.ndarray)
        assert omv.initvec.size == 0
        assert omv.initvec.dtype == np.float32

    # Initializing with skip_calc=True prevents recalculation in get_mean_stdev
    def test_skip_calc_prevents_recalculation(self):
        import numpy as np
        from src.OnlineMeanVar import OnlineMeanVar
    
        # Create instance with skip_calc=True
        omv = OnlineMeanVar(skip_calc=True)
    
        # Set some values manually
        omv.mean = 10.0
        omv.std = 2.0
    
        # Add some data
        omv.append(np.array([5.0, 15.0, 10.0], dtype=np.float32))
    
        # Get mean and std - should return the manually set values without recalculation
        mean, std = omv.get_mean_stdev()
    
        # Verify that the values weren't recalculated
        assert mean == 10.0
        assert std == 2.0
    
        # Verify that _flush wasn't called (buffer should still contain our values)
        assert omv.initvec.size == 3

    # Appending enough values to exceed INITLEN triggers automatic processing
    def test_append_triggers_automatic_processing(self):
        import numpy as np
        from src.OnlineMeanVar import OnlineMeanVar
    
        # Create instance with default parameters
        omv = OnlineMeanVar(initlen=5)  # Set INITLEN to 5 for testing
    
        # Append values less than INITLEN
        omv.append(np.array([1.0, 2.0, 3.0], dtype=np.float32))
        assert omv.n == 0  # No processing should have occurred
    
        # Append more values to exceed INITLEN
        omv.append(np.array([4.0, 5.0], dtype=np.float32))
    
        # Check if processing occurred
        assert omv.n == 5  # Processing should have occurred
        assert omv.initvec.size == 0  # Buffer should be reset

    # Appending values to an empty buffer stores them for later processing
    def test_append_to_empty_buffer(self):
        import numpy as np
        from src.OnlineMeanVar import OnlineMeanVar
    
        # Create instance with default parameters
        omv = OnlineMeanVar()
    
        # Append values to the empty buffer
        values_to_append = np.array([1.0, 2.0, 3.0], dtype=np.float32)
        omv.append(values_to_append)
    
        # Check if values are stored in the buffer
        assert np.array_equal(omv.initvec, values_to_append)

    # Calling get_mean_stdev after appending values returns correct mean and standard deviation
    def test_get_mean_stdev_after_appending_values(self):
        import numpy as np
        from src.OnlineMeanVar import OnlineMeanVar
    
        # Create instance with default parameters
        omv = OnlineMeanVar()
    
        # Append values
        values = np.array([1.0, 2.0, 3.0, 4.0, 5.0], dtype=np.float32)
        omv.append(values)
    
        # Get mean and standard deviation
        mean, std = omv.get_mean_stdev()
    
        # Calculate expected mean and standard deviation
        expected_mean = np.mean(values)
        expected_std = np.std(values, ddof=1)
    
        # Assert the results
        assert np.isclose(mean, expected_mean, atol=1e-6)
        assert np.isclose(std, expected_std, atol=1e-6)

    # Multiple append operations accumulate statistics correctly over time
    def test_multiple_appends_accumulate_statistics(self):
        import numpy as np
        from src.OnlineMeanVar import OnlineMeanVar
    
        # Create instance with default parameters
        omv = OnlineMeanVar()
    
        # Append first batch of values
        omv.append(np.array([1.0, 2.0, 3.0]))
        mean, std = omv.get_mean_stdev()
    
        # Check statistics after first append
        assert omv.n == 3
        assert np.isclose(mean, 2.0)
        assert np.isclose(std, 1.0)
    
        # Append second batch of values
        omv.append(np.array([4.0, 5.0, 6.0]))
        mean, std = omv.get_mean_stdev()
    
        # Check statistics after second append
        assert omv.n == 6
        assert np.isclose(mean, 3.5)
        assert np.isclose(std, 1.8708, atol=1e-4)