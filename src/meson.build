utils_src = [ 'utils.cpp', 'utils.hpp' ]

build_index_src = [ 'goldrush_edit_index.cpp' ] + utils_src
build_targeted_bfs_src = [ 'goldrush_edit_targeted_bfs.cpp' ] + utils_src

executable('goldrush-edit-index',
            build_index_src,
            dependencies : deps,
            install : true)
executable('goldrush-edit-targeted-bfs',
            build_targeted_bfs_src,
            dependencies: deps,
            install : true)