
from setuptools import setup, find_packages

import re
import ast


# version parsing from __init__ pulled from Flask's setup.py
# https://github.com/mitsuhiko/flask/blob/master/setup.py
_version_re = re.compile(r'__version__\s+=\s+(.*)')

with open('q2_dbbact/__init__.py', 'rb') as f:
    hit = _version_re.search(f.read().decode('utf-8')).group(1)
    version = str(ast.literal_eval(hit))

classifiers = [
    'Development Status :: 4 - Beta',
    'License :: OSI Approved :: BSD License',
    'Environment :: Plugins',
    'Topic :: Software Development :: Libraries :: Application Frameworks',
    'Topic :: Scientific/Engineering',
    'Topic :: Scientific/Engineering :: Bio-Informatics',
    'Programming Language :: Python',
    'Programming Language :: Python :: 3',
    'Programming Language :: Python :: 3.6',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
    'Operating System :: Unix',
    'Operating System :: POSIX',
    'Operating System :: MacOS :: MacOS X',
    'Operating System :: Microsoft :: Windows']


description = 'A qiime2 (https://qiime2.org/) plugin for dbBact (http://dbbact.org) annotations of microbiome experiments'

with open('README.md') as f:
    long_description = f.read()

keywords = 'microbiome qiime2 dbbact database analysis bioinformatics',

setup(
    name="q2-dbbact",
    version=version,
    packages=find_packages(),
    license='BSD',
    long_description=long_description,
    long_description_content_type="text/markdown",
    keywords=keywords,
    classifiers=classifiers,
    author="dbBact team",
    author_email="info@dbbact.org",
    url="http://dbbact.org",
    description=description,
    python_requires='>=3.6',
    entry_points={
        "qiime2.plugins":
        ["q2-dbbact=q2_dbbact.plugin_setup:plugin"]
    },
    install_requires=[
        'wordcloud',
        'calour',
        'dbbact-calour'
    ]
)
