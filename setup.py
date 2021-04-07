import os
from setuptools import setup, find_packages
from src.version import __version__

cur_path = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(cur_path, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(name='SVision',
      version=__version__,

      description='SV/CSV callers',
      long_description=long_description,

      url='',
      author='Jiadong Lin, Songbo Wang',
      author_email='jiadong66@stu.xjtu.edu.cn, songbowang125@163.com',

      license='GPLv3',
      classifiers=[
      'Operating System :: POSIX :: Linux',
      'Topic :: Scientific/Engineering :: Bio-Informatics',
      'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
      'Programming Language :: Python :: 3.6'
      ],
      keywords=['SVision', 'deep learning', 'sv', 'csv', 'long read'],

      packages = ['src', 'src/collection', 'src/network', 'src/segmentplot'],
      data_files = [("", ["LICENSE"])],

      zip_safe=False,
      python_requires='>=3.6',
      # install_requires=['scipy>=1.3.0', 'pysam>=0.15.3', 'tensorflow>=1.14.0', 'opencv-python>=4.1.1'],
      scripts=['SVision.py'],

      entry_points={'console_scripts': ['SVision=SVision:main']}

      )
