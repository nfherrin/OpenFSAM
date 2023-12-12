project(
  'OpenFSAM',
  'fortran',
  license: 'MIT',
  license_files: 'LICENSE',
  default_options: [
    'buildtype=debugoptimized',
    'default_library=both',
  ],
)
install = not (meson.is_subproject() and get_option('default_library') == 'static')

# Grab OpenFSAM Sources
srcs = []
subdir('src')

# Build/Install OpenFSAM
openfsam_lib = library(
  meson.project_name(),
  sources: srcs,
  install: install,
)

pkg = import('pkgconfig')
pkg.generate(
  openfsam_lib,
  description: 'An open source Fortran based simulated annealing utility',
)

# Export OpenFSAM as a Dependency for Other Projects (and Examples)
openfsam_inc = openfsam_lib.private_dir_include()
openfsam_dep = declare_dependency(
  link_with: openfsam_lib,
  include_directories: openfsam_inc,
)

# Build the Examples
subdir('examples')
