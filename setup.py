from setuptools import setup

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='ldsc',
      version='3.0.2',
      description='LD Score Regression (LDSC)',
      long_description=long_description,
      long_description_content_type='text/markdown',
      url='http://github.com/bulik/ldsc',
      author='Brendan Bulik-Sullivan and Hilary Finucane',
      author_email='',
      license='GPLv3',
      packages=['ldscore'],
      scripts=['ldsc.py', 'munge_sumstats.py', 'make_annot.py'],
      install_requires = [
            'bitarray==2.6.0',
            'nose==1.3.7',
            'numpy==1.23.3',
            'pandas==1.5.0',
            'pybedtools==0.9.1',
            'pysam==0.19.1',
            'python-dateutil==2.8.2',
            'pytz==2022.4',
            'scipy==1.9.2',
            'six==1.16.0'
      ]
)
