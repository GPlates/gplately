# To install locally: python setup.py build && python setup.py install
# (If there are problems with installation of the documentation, the
#  egg file may be out of sync and will need to be manually deleted 
#  - see error message for details of the corrupted zip file. )
#
# To push a version through to pip.
#  - Make sure it installs correctly locally as above
#  - Update the version information in this file
#  - python setup.py sdist upload -r pypitest  # for the test version
#  - python setup.py sdist upload -r pypi      # for the real version
# With twine:
#  - python setup.py sdist
#  - twine upload dist/*
#
# (see http://peterdowns.com/posts/first-time-with-pypi.html)


from setuptools import dist, setup, find_packages

# try: 
#     from distutils.command import bdist_conda
# except ImportError:
#     pass

from os import path
import io
import os
import subprocess
import platform 

link_args = []
 
if "Windows" in platform.system():
    link_args = ["-static"]

# in development set version to none and ...
PYPI_VERSION = "1.0"  # Note: don't add any dashes if you want to use conda, use b1 not .b1 

# Return the git revision as a string (from numpy)

def git_version():
    
    def _minimal_ext_cmd(cmd):
        # construct minimal environment
        env = {}
        for k in ['SYSTEMROOT', 'PATH']:
            v = os.environ.get(k)
            if v is not None:
                env[k] = v
        # LANGUAGE is used on win32
        env['LANGUAGE'] = 'C'
        env['LANG'] = 'C'
        env['LC_ALL'] = 'C'
        out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=env).communicate()[0]
        return out

    try:
        out = _minimal_ext_cmd(['git', 'rev-parse', '--short', 'HEAD'])
        GIT_REVISION = out.strip().decode('ascii')
    except OSError:
        GIT_REVISION = "Unknown"

    return GIT_REVISION


if PYPI_VERSION is None:
    PYPI_VERSION = git_version()


this_directory = path.abspath(path.dirname(__file__))
with io.open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


if __name__ == "__main__":
    setup(name = 'gplately',
          author            = "Ben Mather",
          author_email      = "ben.mather@sydney.edu.au",
          url               = "https://github.com/GPlates/gplately",
          version           = PYPI_VERSION,
          description       = "Object-orientated Python interface to pyGPlates for plate tectonic reconstructions",
          long_description  = long_description,
          long_description_content_type='text/markdown',
          install_requires  = ['numpy>=1.16.0',
                               'scipy>=1.0.0',
                               'shapely',
                               'matplotlib',
                               'cartopy',
                               'PlateTectonicTools',
                               'pooch',
                               'tqdm',
                               'netcdf4',
                               'rasterio',
                               'geopandas',
                               'stripy',
                               ],
          packages          = ['gplately'],
          package_data      = {'gplately': ['Notebooks/*ipynb', # Worked Examples is not currently used
                                            'Notebooks/Data/*'] },
          include_package_data = True,
          classifiers       = ['Programming Language :: Python :: 3',
                               'Programming Language :: Python :: 3.5',
                               'Programming Language :: Python :: 3.6',
                               'Programming Language :: Python :: 3.7',
                               'Programming Language :: Python :: 3.8',
                               'Programming Language :: Python :: 3.9',
                               'Programming Language :: Python :: 3.10',
                               'Programming Language :: Python :: 3.11'
                               ]
          )
