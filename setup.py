import os
from setuptools import setup, find_packages
from src.version import __version__

cur_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(cur_path, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
      name='SVision',
      version=__version__,

      description='SV/CSV callers',
      long_description=long_description,

      url='https://github.com/xjtu-omics/SVision',
      # download_url = 'https://github.com/xjtu-omics/SVision/archive/refs/tags/v1.2.2-beta.tar.gz',

      author='Jiadong Lin, Songbo Wang',
      author_email='jiadong66@stu.xjtu.edu.cn, songbowang125@163.com',

      license='GPLv3',
      classifiers=[
      'Operating System :: POSIX :: Linux',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      'Programming Language :: Python :: 3.6'
      ],
      keywords=['SVision', 'Deep learning', 'Complex structural variants', 'Structural variants', 'Single moleculo sequencing'],

      packages = ['src', 'src/collection', 'src/network', 'src/segmentplot'],
      data_files = [("", ["LICENSE"])],

      zip_safe=False,
      python_requires='>=3.6',
      install_requires=['scipy', 'beautifulsoup4', 'numpy==1.16'],
      scripts=['SVision.py'],

      entry_points={'console_scripts': ['SVision=SVision:main']}

      )
