from setuptools import setup, find_packages

requirements = [
      'python',
      'h5py',
      'biopython',
      'mafft',
      'matplotlib',
      'numpy',
      'scipy',
      'nanopolish',
      'minimap2',
      'pandas',
      'seaborn',
      'psutil',
      'hdf5plugin',
      'ont_vbz_hdf_plugin'
      ]

setup(name='Magnipore',
      version='1.0.0',
      description='Compares two ONT sequences samples for differential signals cause by mutations and modifications',
      author='Jannes Spangenberg',
      author_email='jannes.spangenberg@uni-jena.de',
      url='https://github.com/JannesSP/magnipore',
      license='GNU General Public License v3.0',
      packages=find_packages(),
      py_modules=['magnipore'],
      package_data={
        'magnipore':['README.md','LICENSE'],
      },
      entry_points={
            'console_scripts':[
                  'magnipore=magnipore.magnipore:main',
                  'nanosherlock=magnipore.nanosherlock:main',
            ]
      },
      install_requires=requirements,
      keywords=['Magnipore', 'magnipore', 'nanosherlock', 'ONT', 'Oxford Nanopore Technologies', 'MinION', 'Direct RNA Sequencing'],
     )
