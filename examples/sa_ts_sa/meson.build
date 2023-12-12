sa_ts_sa_srcs = files(
  'globals.f90',
  'infuncs.f90',
  'outfuncs.f90',
  'sa_ts_sa.f90',
  'travel_sales.f90',
)

sa_ts_sa_exec = executable(
  'sa_ts_sa',
  'main.f90',
  sources: sa_ts_sa_srcs,
  dependencies: openfsam_dep,
)
