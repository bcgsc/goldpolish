src_files = [ 'FastaIndex.h', 'FastaInterleave.h', 'FastaReader.cpp', 'FastaReader.h', 'FastaWriter.cpp', 'FastaWriter.h', 'FastaConcat.h', 'Options.h' ]

datalayer_lib = static_library('datalayer',
        src_files,
        dependencies : deps + [ common_dep ],
        install : false, include_directories : [ '..' ])
datalayer_dep = declare_dependency(link_with : datalayer_lib)