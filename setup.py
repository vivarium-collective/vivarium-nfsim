import setuptools
from distutils.core import setup

with open("README.md", 'r') as readme:
    long_description = readme.read()

setup(
    name='vivarium-nfsim',
    version='0.0.1',
    packages=[
        'vivarium_nfsim',
        'vivarium_nfsim.processes',
        'vivarium_nfsim.experiments',
        'vivarium_nfsim.models',
    ],
    author='Eran Agmon',
    author_email='agmon.eran@gmail.com',
    url='https://github.com/vivarium-collective/vivarium-nfsim',
    license='Apache-2.0',
    short_description='A process-bigraph wrapper for BioNetGen/NFSim.',
    long_description=long_description,
    long_description_content_type='text/markdown',
    package_data={},
    include_package_data=True,
    install_requires=[
        'process-bigraph',
        'bionetgen',
        'pytest',
    ],
)
