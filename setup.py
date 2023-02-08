from setuptools import setup, find_packages
import setuptools
from setuptools.extension import Extension
from Cython.Build import cythonize
import numpy
from distutils import ccompiler


# This was stolen from pybind11
# https://github.com/pybind/python_example/blob/master/setup.py
# As of Python 3.6, CCompiler has a `has_flag` method.
# cf http://bugs.python.org/issue26689


def has_flag(compiler, flagname):
    """Return a boolean indicating whether a flag name is supported on
    the specified compiler.
    """
    import tempfile

    with tempfile.NamedTemporaryFile('w', suffix='.cpp') as f:
        f.write('int main (int argc, char **argv) { return 0; }')
        try:
            compiler.compile([f.name], extra_postargs=[flagname])
        except setuptools.distutils.errors.CompileError:
            return False
    return True


def cpp_flag(compiler, flags):
    """Return the -std=c++[11/14/17] compiler flag.
    The newer version is prefered over c++11 (when it is available).
    """
    for flag in flags:
        if has_flag(compiler, flag):
            return flag


def get_extra_args():
    compiler = ccompiler.new_compiler()
    extra_compile_args = []

    flags = ['-std=c++20', '-std=c++17', '-std=c++14', '-std=c++11']

    f = cpp_flag(compiler, flags)
    if not f:
        raise RuntimeError("Invalid compiler")
    extra_compile_args.append(f)

    flags = ['--stdlib=libc++']
    f = cpp_flag(compiler, flags)
    if f:
        extra_compile_args.append(f)
    # flags = ['-W#warnings']
    # f = cpp_flag(compiler, flags)
    # if f:
    #     extra_compile_args.append(f)
    return extra_compile_args


extras = get_extra_args()
print("Extra compiler args ", extras)


ext_modules = []
for item in ["io_funcs", "input_stream_alignments", "pairing", "samclips"]:

    ext_modules.append(Extension(f"dodi.{item}",
                                 [f"dodi/{item}.pyx"],
                                 library_dirs=[numpy.get_include()],
                                 extra_compile_args=extras,
                                 language="c++"))


print("Found packages", find_packages(where="."))
setup(
    name="dodi",
    version='0.4.3',
    python_requires='>=3.7',
    install_requires=[
            'cython',
            'click',
            'numpy',
            'ncls'
        ],
    packages=find_packages(where="."),
    ext_modules=cythonize(ext_modules),
    include_dirs=[numpy.get_include()],
    include_package_data=True,
    entry_points='''
        [console_scripts]
        dodi=dodi.main:dodi_aligner
    ''',
)
