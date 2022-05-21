src_files = [ 'alignGlobal.cc', 'alignGlobal.h', 'dialign.cpp', 'dialign.h','dna_diag_prob.cc', 'smith_waterman.cpp', 'smith_waterman.h', 'Options.h' ]

align_lib = static_library('align', src_files, dependencies : deps + [ common_dep, dialign_dep ], install : false, include_directories : [ '..', '../Common' ])
align_dep = declare_dependency(link_with : align_lib)