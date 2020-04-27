from setuptools import setup, find_packages
from os import path, environ

cur_dir = path.abspath(path.dirname(__file__))

with open(path.join(cur_dir, 'requirements.txt'), 'r') as f:
    requirements = f.read().split()

setup(
    name='pytao',
    version = '0.1.0',
    packages=find_packages(),  
    package_dir={'pytao':'pytao'},
    url='https://www.classe.cornell.edu/bmad/tao.html',
    long_description=open('README').read(),
    long_description_content_type='text/markdown',
    install_requires=requirements,
    include_package_data=True,
    python_requires='>=3.6'
)
