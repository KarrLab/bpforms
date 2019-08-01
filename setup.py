import setuptools
try:
    import pkg_utils
except ImportError:
    import pip._internal
    pip._internal.main(['install', 'pkg_utils'])
    import pkg_utils
import os

name = 'bpforms'
dirname = os.path.dirname(__file__)
package_data = {
    name: [
        'VERSION',
        'grammar.lark',
        'alphabet/*.yml',  # alphabets
        'config/*.cfg',
        'web/*',  # website
        'web/css/*',  # website
        'web/js/*',  # website
    ],
}

# get package metadata
md = pkg_utils.get_package_metadata(dirname, name, package_data_filename_patterns=package_data)

# install package
setuptools.setup(
    name=name,
    version=md.version,
    description='Unambiguous representation of modified DNA, RNA, and proteins',
    long_description=md.long_description,
    url="https://www.bpforms.org",
    download_url='https://github.com/KarrLab/' + name,
    author="Karr Lab",
    author_email="info@karrlab.org",
    license="MIT",
    keywords='DNA,RNA,protein,post-transcriptional modification,post-translational modification,proteoform,phosphorylation,methylation',
    packages=setuptools.find_packages(exclude=['tests', 'tests.*']),
    package_data=md.package_data,
    install_requires=md.install_requires,
    extras_require=md.extras_require,
    tests_require=md.tests_require,
    dependency_links=md.dependency_links,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
        'Natural Language :: English',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
    ],
    entry_points={
        'console_scripts': [
            'bpforms = bpforms.__main__:main',
        ],
    },
)
