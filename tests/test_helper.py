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
