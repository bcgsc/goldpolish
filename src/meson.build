common = [ 'utils.cpp', 'utils.hpp', 'seqindex.cpp', 'seqindex.hpp', 'mappings.cpp', 'mappings.hpp', 'fn_name.hpp' ]

build_index_src = [ 'goldrush_edit_index.cpp' ] + common
build_targeted_bfs_src = [ 'goldrush_edit_targeted_bfs.cpp' ] + common
build_hold_src = [ 'goldrush_edit_hold.cpp' ] + common

executable('goldrush-edit-index',
            build_index_src,
            dependencies : deps,
            install : true)
executable('goldrush-edit-targeted-bfs',
            build_targeted_bfs_src,
            dependencies: deps,
            install : true)
executable('goldrush-edit-hold',
            build_hold_src,
            dependencies: deps,
            install : true)