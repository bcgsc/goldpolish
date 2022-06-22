project('sealer-gr', [ 'cpp', 'c' ],
        default_options : [ 'cpp_std=c++17', 'warning_level=0', 'werror=false' ])

compiler = meson.get_compiler('cpp')

threads_dep = dependency('threads')
openmp_dep = dependency('openmp')
btllib_dep = compiler.find_library('btllib')
boost_deps = [ compiler.find_library('boost_graph') ]

deps = [ threads_dep, openmp_dep, btllib_dep ] + boost_deps

subdir('Common')
subdir('dialign')
subdir('Align')
subdir('DataLayer')
subdir('Assembly')

executable('sealer-gr',
           'sealer.cpp',
            dependencies : deps + [ common_dep, align_dep, datalayer_dep ],
            install : true)