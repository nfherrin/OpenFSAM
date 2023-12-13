traveling_sales_general_srcs = files(
  'globals.f90',
  'infuncs.f90',
  'outfuncs.f90',
  'travel_sales.f90',
)

traveling_sales_general_exec = executable(
  'traveling_sales_general',
  'main.f90',
  sources: traveling_sales_general_srcs,
  dependencies: openfsam_dep,
)
