"""
Setup Xepto 50
"""
import sys
try:
    from setuptools import setup, find_packages
except ImportError:
    raise ImportError("Please install `setuptools`")

if not (sys.version_info[0] == 3 and sys.version_info[1] >= 9):
    raise RuntimeError(
                'alphaML requires Python 3.9 or higher.\n'
                'You are using Python {0}.{1}'.format(
                    sys.version_info[0], sys.version_info[1])
                )

# main setup command
setup(
    name='xepto50',
    version='0.0.2', # Major.Minor.Patch
    author='Julhash Kazi',
    author_email='xepto50@kazilab.se',
    url='https://www.kazilab.se',
    description='Calculate drug sensitivity in batch mode!',
    license='Apache-2.0',
    install_requires=[
        'IPython',
        'matplotlib',
        'numpy',
        'scipy',
        'pandas',
        'scikit-learn',
        'lmfit'
    ],
    extras_require={
        'gui': ['PyQt5']
    },
    entry_points={
        'console_scripts': [
            'xepto50=xepto50:main',  # To maps "xepto50" command to the main function
        ],
    },
    platforms='any',
    packages=find_packages()
)
