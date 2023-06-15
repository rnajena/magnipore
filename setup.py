from setuptools import setup
import versioneer

requirements = [
      'h5py>=3.7',
      'biopython>=1.80',
      'mafft>=7.508',
      'matplotlib>=3.6',
      'numpy>=1.23',
      'scipy>=1.9',
      'nanopolish>=0.14',
      'minimap2>=2.24',
      'samtools>=1.0',
      'pandas>=1.5',
      'seaborn>=0.12',
      'psutil>=5.9',
      'hdf5>=1.12',
      'hdf5plugin>=3.3',
      'ont_vbz_hdf_plugin>=1.0',
      'pytest>=7.1',
      ]

setup(
    name='magnipore',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Compares two ONT sequences samples for differential signals cause by mutations and modifications.",
    license="GNUv3",
    author="Jannes Spangenberg",
    author_email='jannes.spangenberg@uni-jena.de',
    url='https://github.com/JannesSP/magnipore',
    packages=['magnipore'],
    entry_points={
        'console_scripts':[
            'magnipore=magnipore.magnipore:main',
            'nanosherlock=magnipore.nanosherlock:main',
        ]
    },
    install_requires=requirements,
    keywords=['Magnipore', 'magnipore', 'nanosherlock', 'ONT', 'Oxford Nanopore Technologies', 'MinION', 'Direct RNA Sequencing'],
    classifiers=[
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ]
)
