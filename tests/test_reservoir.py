import numpy as np

class TestReservoir:
    
    # Adding elements to an empty reservoir until it's full
    def test_adding_elements_until_full(self):
        
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
    #     
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