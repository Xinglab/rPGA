import sys, os

# should be able to safely do this now.
from setuptools import setup, find_packages

setup(name='rPGA',
      version='2.0.0',
      packages = find_packages('src'),  # include all packages under src
			package_dir = {'':'src'},   # all distutils packages are under src
      entry_points={'console_scripts': ['rPGA=rPGA.scripts.rPGA:main']},
			description = 'rPGA',
			author = 'Shayna R. Stein, Emad Bahrami-Samani',
			author_email = 'sstein93@ucla.edu',
			url = 'https://github.com/xinglab/rPGA',
			download_url = 'https://github.com/xinglab/rPGA/tarball/rPGA_2.0.0',
			license='GPL3',
			keywords = [],
			classifiers = [],
)
