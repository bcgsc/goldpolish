common = [ 'utils.cpp', 'utils.hpp', 'seqindex.cpp', 'seqindex.hpp', 'mappings.cpp', 'mappings.hpp', 'fn_name.hpp' ]

build_index_src = [ 'goldpolish_index.cpp' ] + common
build_targeted_bfs_src = [ 'goldpolish_targeted_bfs.cpp' ] + common
build_hold_src = [ 'goldpolish_hold.cpp' ] + common

executable('goldpolish-index',
            build_index_src,
            dependencies : deps,
            install : true)
executable('goldpolish-targeted-bfs',
            build_targeted_bfs_src,
            dependencies: deps,
            install : true)
executable('goldpolish-hold',
            build_hold_src,
            dependencies: deps,
            install : true)
