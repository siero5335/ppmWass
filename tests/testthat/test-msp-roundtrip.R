test_that("write_msp -> read_msp preserves entry IDs and non-empty spectra", {
  df_spec <- tibble::tibble(
    id = c("cmpd_A", "cmpd_B"),
    RI = c(1000, 1100),
    known = c("T", "T")
  )

  frag_list <- list(
    cmpd_A = rbind(c(50.0000, 100), c(73.0474, 250), c(91.0548, 500)),
    cmpd_B = rbind(c(43.0548, 300), c(57.0704, 120), c(147.0654, 200))
  )

  tmp <- tempfile(fileext = ".msp")
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)

  ppmWass::write_msp(df_spec, frag_list, tmp)
  lib <- ppmWass::read_msp(tmp, progress = 0)

  expect_true(all(c("cmpd_A", "cmpd_B") %in% lib$name))
  expect_true(all(nchar(lib$spectrum) > 0))

  df_int <- ppmWass::convert_msp_to_internal(lib)
  expect_equal(sort(df_int$id), sort(c("cmpd_A", "cmpd_B")))
  expect_true(all(nchar(df_int$ei) > 0))
})

test_that("read_msp splits adjacent records when peak blocks lack blank separators", {
  tmp <- tempfile(fileext = ".msp")
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)

  writeLines(c(
    "Name: cmpd_A",
    "Comment: no Num Peaks field",
    "50 100",
    "73 200",
    "Name: cmpd_B",
    "Comment: second record starts immediately",
    "60 120",
    "91 180"
  ), tmp)

  lib <- ppmWass::read_msp(tmp, progress = 0)
  expect_equal(lib$name, c("cmpd_A", "cmpd_B"))
  expect_equal(lib$spectrum, c("50:100 73:200", "60:120 91:180"))
})

test_that("build_spectra_from_msp processes spectra through shared pipeline", {
  df_spec <- tibble::tibble(
    id = c("cmpd_A", "cmpd_B"),
    RI = c(1000, 1100),
    known = c("T", "T")
  )

  frag_list <- list(
    cmpd_A = rbind(c(50.0000, 100), c(73.0474, 250), c(91.0548, 500)),
    cmpd_B = rbind(c(43.0548, 300), c(57.0704, 120), c(147.0654, 200))
  )

  tmp <- tempfile(fileext = ".msp")
  on.exit({
    if (file.exists(tmp)) unlink(tmp)
  }, add = TRUE)

  ppmWass::write_msp(df_spec, frag_list, tmp)
  res <- ppmWass::build_spectra_from_msp(
    tmp,
    params = ppmWass::eihrms_default_params(),
    require_ri = FALSE,
    progress = FALSE
  )

  expect_equal(sort(names(res$frag_list)), sort(df_spec$id))
  expect_equal(sort(names(res$loss_list)), sort(df_spec$id))
  expect_true(all(res$df_spec$known == "T"))
})
