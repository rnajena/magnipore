class TestCli:
        
    # Function prints the correct version string with the expected format
    def test_prints_correct_version_format(self, mocker, capsys):
        # Arrange
        mock_version = "v1.2.3"
        mocker.patch("src.cli.__version_str__", mock_version)
        mock_exit = mocker.patch("sys.exit")
    
        # Act
        from src.cli import print_version
        print_version()
    
        # Assert
        captured = capsys.readouterr()
        assert captured.out == f"magnipore {mock_version}\n"
        mock_exit.assert_called_once_with(0)

    # When subtool is not in script_mapping and not a help/version flag, error message is shown
    def test_invalid_subtool_shows_error(self, mocker, capsys):
        import sys
        # Arrange
        mocker.patch('sys.argv', ['magnipore', 'invalid_subtool'])
        mocker.patch('sys.exit')
    
        # Act
        from src.cli import main
        main()
    
        # Assert
        captured = capsys.readouterr()
        assert "Error: 'invalid_subtool' not found" in captured.out
        sys.exit.assert_called_once_with(404)