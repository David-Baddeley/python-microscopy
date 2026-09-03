
install_dir = py.get_install_dir() / 'PYME/simulation'

py_sources = files(
    '__init__.py',
    'wormlike.py',
    'wormlike2.py',
    'wormlike3.py',
    'pointsets.py',
    'wormlike_curled.py',
    'locify.py',
    'image_gen.py',
    'cameranoise.py',
    'fluorophores.py',
)

py.install_sources(py_sources, subdir:'PYME/simulation')

data_files = files(
    'illuminate.pyx',
)

install_data(data_files, install_dir: install_dir)

py.extension_module(
    'illuminate',
    [
        'illuminate.pyx',
    ],
    include_directories: [np_include_dir],
    subdir: 'PYME/simulation',
    dependencies: py.dependency(),
    c_args: ['-O3', '-fno-exceptions', '-ffast-math'],
    install: true,
    #cpp_args: ['-std=c++17'],
)
