from setuptools import find_packages, setup
setup(
    name='slcode',
    packages=find_packages(include=['SL_C0de']),
    version='0.3.5',
    description='A library wich solve the sea level equation for a radial symetrical spherical earth',
    author='HENRY Adrien',
    license='MIT',
    install_requires=['stripy','numpy','pyshtools','scipy','matplotlib','plotly'],
    setup_requires=['pytest-runner'],
    tests_require=['pytest'],
    test_suite='tests',
)