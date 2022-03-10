from setuptools import setup, find_packages

# file containing package dependencies names and versions
with open('requirements.txt') as f:
    required = f.read().splitlines()
print(required)
setup(name='GeneQuesst',
      version='1.0',
      description='Genome assembly tool optomized for gene searching',
      author='Erik Serrano',
      author_email='erik.serrano@cuanschutz.edu',
      packages=find_packages(),
      install_requires=required
     )