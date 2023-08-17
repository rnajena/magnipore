from setuptools import setup
import versioneer

# conda dependencies
requirements = [
    'h5py>=3.7',
    'biopython>=1.80',
    'mafft>=7.508',
    'matplotlib>=3.6',
    'numpy>=1.21',
    'scipy>=1.9',
    'minimap2>=2.24',
    'samtools>=1.0',
    'pandas>=1.5',
    'seaborn>=0.12',
    'psutil>=5.9',
    'hdf5>=1.12',
    'hdf5plugin>=3.3',
    'ont_vbz_hdf_plugin>=1.0',
    'pytest>=7.1',
    'gzip>=1.12',
    'f5c>=1.2',
    'read5>=1.2'
]

setup(
    name='magnipore',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="Compares two ONT sequences samples for differential signals cause by mutations and modifications.",
    license="GNU General Public License v3",
    author="Jannes Spangenberg",
    author_email='jannes.spangenberg@uni-jena.de',
    url='https://github.com/JannesSP/magnipore',
    packages=['magnipore'],
    entry_points={
        'console_scripts':[
            'magnipore=magnipore.magnipore:main',
            'nanosherlock=magnipore.nanosherlock:main',
            'magnifilter=magnipore.magnifilter:main',
            'magniplot=magnipore.magniplot:main',
        ]
    },
    install_requires=requirements,
    python_requires='>=3.8,<3.11',
    keywords=['Magnipore', 'magnipore', 'nanosherlock', 'ONT', 'Oxford Nanopore Technologies', 'MinION', 'Direct RNA Sequencing'],
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Natural Language :: English',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ]
)
