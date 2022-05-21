src_files = [ 'AssemblyAlgorithms.cc', 'AssemblyAlgorithms.h', 'BranchGroup.h', 'BranchRecord.h', 'BranchRecordBase.h', 'DBG.h', 'DotWriter.h', 'Options.cc', 'Options.h', 'SequenceCollection.h', 'VertexData.h', 'AdjacencyAlgorithm.h', 'AssembleAlgorithm.h', 'BubbleAlgorithm.h', 'CoverageAlgorithm.h', 'ErodeAlgorithm.h', 'LoadAlgorithm.h', 'SeqExt.h', 'SplitAlgorithm.h', 'TrimAlgorithm.h' ]

assembly_lib = static_library('assembly',
        src_files,
        dependencies : deps + [ common_dep, datalayer_dep ],
        install : false,
        include_directories : '..')
assembly_dep = declare_dependency(link_with : assembly_lib)