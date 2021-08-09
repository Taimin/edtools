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
    version='1.0.1',
    description='Collection of tools for automated processing and clustering of electron diffraction data.',
    python_requires='>=3.6.1',
    project_urls={
        'documentation': 'http://github.com/instamatic-dev/edtools',
        'homepage': 'http://github.com/instamatic-dev/edtools',
        'repository': 'http://github.com/instamatic-dev/edtools'},
    author='Stef Smeets',
    author_email='s.smeets@esciencecenter.nl',
    license='GPL-3.0-only',
    keywords='electron-diffraction microed xds pipeline cluster-analysis',
    classifiers=[
            'Development Status :: 5 - Stable',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
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
            'edtools.group_reflections = edtools.group_reflections:main',
            'edtools.find_rotation_axis = edtools.find_rotation_axis:main',
            'edtools.find_beam_center  = edtools.find_beam_center:main',
            'edtools.scattering_factor  = edtools.scattering_factor:main']},
    packages=['edtools'],
    package_dir={
        '': '.'},
    package_data={
        'edtools': ['*.yaml']},
    install_requires=[
        'matplotlib==3.*,>=3.2.1',
        'numpy==1.*,>=1.18.2',
        'pandas==1.*,>=1.0.3',
        'scipy==1.*,>=1.4.1',
        'uncertainties==3.*,>=3.1.2',
        'lmfit>=1.0.0'],
    extras_require={
        'dev': [
            'check-manifest',
            'pre-commit']},
)
