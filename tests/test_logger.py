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
