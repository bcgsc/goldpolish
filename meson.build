project('goldrush-edit', 'cpp',
        version : '0.1.0',
        license : 'GPL3',
        default_options : [ 'cpp_std=c++17', 'warning_level=3', 'werror=true' ])

# Compiler and global flags
# ===========================================================

compiler = meson.get_compiler('cpp')
compiler_id = compiler.get_id()
compiler_version = compiler.version()

if compiler_id == 'gcc' and compiler_version.version_compare('<7.0.0')
  error('GCC ' + compiler_version + ' is unsupported.')
endif
if compiler_id == 'clang' and compiler_version.version_compare('<10.0.0')
  error('Clang ' + compiler_version + ' is unsupported.')
endif

# Build dependencies and subprojects
# ===========================================================

threads_dep = dependency('threads')
openmp_dep = dependency('openmp')
btllib_dep = compiler.find_library('btllib')

deps = [ threads_dep, openmp_dep, btllib_dep ]

# Source files, scripts
# ===========================================================

subdir('src')
subdir('scripts')