import os
import sys
import shutil
import platform
import subprocess
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext


class CustomBuildExt(build_ext):
    def run(self):
        ext = self.extensions[0]
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        out_path = os.path.join(extdir, 'al2co')
        os.makedirs(extdir, exist_ok=True)

        compiled = False
        try:
            print(f"Compiling al2co.c into {out_path}...")
            subprocess.check_call([
                'gcc', '-O3', '-fopenmp',
                '-Wno-implicit-int', '-Wno-implicit-function-declaration',
                'al2co.c', '-o', out_path, '-lm'
            ])
            compiled = True
            print("Compilation succeeded.")
        except Exception as e:
            print(f"Warning: Could not compile al2co.c using gcc: {e}")
            print("Falling back to pre-built binary.")

        if not compiled:
            # Copy the pre-built binary that ships with the source tree
            src_binary = os.path.join(os.path.dirname(__file__), 'al2co', 'al2co')
            if os.path.exists(src_binary):
                shutil.copy2(src_binary, out_path)
                print(f"Using pre-built binary from {src_binary}")
            else:
                print(
                    "Warning: No pre-built binary found and compilation failed. "
                    "The al2co executable will be missing."
                )


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
