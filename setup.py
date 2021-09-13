from setuptools import setup, Extension


setup(name='spkmeansmodule',
      version='1.0',
      description='spkmeans calculator',
      ext_modules=[Extension('spkmeansmodule', sources=['spkmeans.c','spkmeansmodule.c'])])