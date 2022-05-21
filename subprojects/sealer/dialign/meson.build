src_files = [ 'alig.c', 'assemble.c', 'diag.c', 'io.c', 'io.h', 'orf.c', 'orf.h', 'parameters.c', 'parameters.h', 'prob.c', 'struct.h', 'translate.c', 'translate.h' ]

dialign_lib = static_library('dialign',
        src_files,
        dependencies : deps,
        install : false, include_directories : [ '..', '../Common' ])
dialign_dep = declare_dependency(link_with : dialign_lib)