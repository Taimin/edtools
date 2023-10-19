# DO NOT EDIT THIS FILE!
# This file has been autogenerated by dephell <3
# https://github.com/dephell/dephell

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


import os.path

readme = ''
here = os.path.abspath(os.path.dirname(__file__))
readme_path = os.path.join(here, 'README.rst')
if os.path.exists(readme_path):
    with open(readme_path, 'rb') as stream:
        readme = stream.read().decode('utf8')


setup(
    long_description=readme,
    name='edtools',
    version='1.0.3',
    description='Collection of tools for automated processing and clustering of electron diffraction data.',
    python_requires='>=3.6.1',
    project_urls={
        'documentation': 'http://github.com/instamatic-dev/edtools',
        'homepage': 'http://github.com/instamatic-dev/edtools',
        'repository': 'http://github.com/instamatic-dev/edtools'},
    author='Stef Smeets',
    author_email='s.smeets@esciencecenter.nl',
    license='BSD-3-clause',
    keywords='electron-diffraction microed xds pipeline cluster-analysis',
    classifiers=[
            'Development Status :: 5 - Stable',
            'License :: OSI Approved :: BSD License',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
    ],
    entry_points={
        'console_scripts': [
            'edtools.autoindex = edtools.autoindex:main',
            'edtools.cluster = edtools.cluster:main',
            'edtools.find_cell = edtools.find_cell:main',
            'edtools.extract_xds_info = edtools.extract_xds_info:main',
            'edtools.run_pointless = edtools.run_pointless:main',
            'edtools.make_xscale = edtools.make_xscale:main',
            'edtools.make_shelx = edtools.make_shelx:main',
            'edtools.update_xds = edtools.update_xds:main',
            'edtools.reflection_tool = edtools.reflection_tool:main',
            'edtools.find_rotation_axis = edtools.find_rotation_axis:main',
            'edtools.find_beam_center  = edtools.find_beam_center:main',
            'edtools.scattering_factor  = edtools.scattering_factor:main',
            'edtools.update_dials  = edtools.update_dials:main',
            'edtools.cif_tools = edtools.cif_tools:main',
            'edtools.autoindex_dials  = edtools.autoindex_dials:main',
            'edtools.extract_dials_info  = edtools.extract_dials_info:main',
            'edtools.dials_to_crystfel  = edtools.dials_to_crystfel:main',]},
    packages=['edtools'],
    package_dir={
        'edtools': 'edtools'},
    package_data={
        'edtools': ['*.yaml',
                    'instrument/*.cif']},
    install_requires=[
        'matplotlib==3.*,>=3.2.1',
        'numpy==1.*,>=1.18.2',
        'openpyxl>=3.0.10',
        'pandas==1.*,>=1.0.3',
        'scipy==1.*,>=1.4.1',
        'uncertainties==3.*,>=3.1.2',
        'lmfit>=1.0.0'],
    extras_require={
        'dev': [
            'check-manifest',
            'pre-commit',
            'bump2version']},
)
