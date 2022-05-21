project('ntedit-gr', 'cpp',
        default_options : [ 'cpp_std=c++17', 'warning_level=3', 'werror=true' ])

compiler = meson.get_compiler('cpp')

threads_dep = dependency('threads')
openmp_dep = dependency('openmp')
btllib_dep = compiler.find_library('btllib')

deps = [ threads_dep, openmp_dep, btllib_dep ]

z_dep = compiler.find_library('z')

executable('ntedit-gr',
            'ntedit.cpp',
            dependencies : deps + z_dep,
            install : true)