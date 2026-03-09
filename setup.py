import os
import sys
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

class CustomBuildExt(build_ext):
    def run(self):
        # Build al2co C binary
        ext = self.extensions[0]
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        out_path = os.path.join(extdir, 'al2co')
        os.makedirs(extdir, exist_ok=True)
        
        try:
            print(f"Compiling al2co.c into {out_path}...")
            # Automatically compiling ensures it drops into the wheel or source tree properly
            subprocess.check_call(['gcc', '-O3', '-fopenmp', 'al2co.c', '-o', out_path, '-lm'])
        except Exception as e:
            print(f"Error compiling al2co.c using gcc: {e}")
            sys.exit(1)

setup(
    name='al2co',
    version='0.1.0',
    description='Conservation score calculation app - wrapped for PyPI',
    author='Antigravity',
    packages=['al2co'],
    ext_modules=[Extension('al2co._dummy', sources=[])],
    cmdclass={'build_ext': CustomBuildExt},
    package_data={'al2co': ['al2co']},
    install_requires=[
        'biopython'
    ]
)
