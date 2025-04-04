class TestMagniplot:

    # Ensure the marker function returns 'D' for uppercase and mixed case mutation contexts.
    def test_marker(self):
        from src.plot import marker
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
        from src.plot import color
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