utils_src = [ 'utils.cpp', 'utils.hpp' ]

build_index_src = [ 'build_index.cpp' ] + utils_src
build_targeted_bfs_src = [ 'build_targeted_bfs.cpp' ] + utils_src

executable('build-index', build_index_src, install : true)
executable('build-targeted-bfs', build_targeted_bfs_src, install : true)