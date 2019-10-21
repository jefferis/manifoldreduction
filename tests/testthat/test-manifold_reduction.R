context("manifold_reduction")

test_that("manifold_reduction works",{
  skip_if_not_installed('nat')
  library(nat)
  x=read.im3d("testdata/testneurons_thresh.nrrd")
  xyz=ind2coord(x)
  xyzm=manifold_reduction(t(xyz), no_iterations = 15)
  expect_known_value(xyzm, file = 'testdata/xyzm.rds')

  # warn for fishy points dimensions
  expect_warning(manifold_reduction(xyz, no_iterations = 1))
})
