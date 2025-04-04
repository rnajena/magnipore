import sys
import subprocess
from os.path import join, dirname, exists

from src.__init__ import __version__, __version_str__

def print_help(exit : int = 0):
    """Display help message for the CLI."""
    help_text = """Usage: magnipore <subtool> [options]

Magnipore: A toolkit for genomic analysis.

Available subtools:
    run       Run the main pipeline
    check     Perform data validation and checks
    filter    Filter sequencing data
    genomic   Analyze genomic features
    plot      Generate plots from data

Options:
    -h, --help    Show this help message and exit
    -v, --version Show the version number and exit
"""
    print(help_text)
    sys.exit(exit)

def print_version():
    """Display version number."""
    print(f"magnipore v{__version__}")
    sys.exit(0)

def main():
    """Main CLI entry point."""
    if len(sys.argv) < 2:
        print_help()

    subtool = sys.argv[1]

    # Handle global options
    if subtool in ["--help", "-h"]:
        print_help()
    elif subtool in ["--version", "-v"]:
        print_version()

    script_mapping = {
        "run": "magnipore.py",
        "check": "check.py",
        "filter": "filter.py",
        "genomic": "genomic.py",
        "plot": "plot.py",
    }

    if subtool in script_mapping:
        script = script_mapping[subtool]
        # Get the absolute path to the script inside the package
        script_path = join(dirname(__file__), script)
        command = [sys.executable, script_path] + sys.argv[2:]

        if not exists(script_path):
            print(f"Error: Script '{script}' not found at {script_path}")
            sys.exit(400)

        sys.exit(subprocess.run(command).returncode)
    else:
        print(f"Error: '{subtool}' not found")
        sys.exit(404)

if __name__ == "__main__":
    main()
