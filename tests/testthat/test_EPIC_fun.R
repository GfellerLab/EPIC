context("Testing EPIC function")

test_that("Test of bad inputs", {
          a <- matrix(1:20, nrow=5)
          expect_error(EPIC(bulk=a), "There are only 0 signature")
          refStr <- "unknownRef"
          expect_error(EPIC(bulk=hoek_data$rawCounts, reference=refStr),
                       paste0(refStr, ".* not part of the allowed references"))
          tRef <- BRef; tRef$sigGenes <- NULL
          expect_error(EPIC(bulk=hoek_data$rawCounts, reference=tRef),
                       "needs to contain .* 'sigGenes'")
          tRef <- BRef;
          tRef$refProfiles.var <- tRef$refProfiles.var[
            nrow(tRef$refProfiles.var):1,]
          expect_error(EPIC(bulk=hoek_data$rawCounts, reference=tRef),
                       "dimnames of .*refProfiles.*refProfiles.var.*same")
})

test_that("Test for correct result on hoek data with default input", {
  testFract <- EPIC(hoek_data$rawCounts)$cellFractions
  expect_equal(testFract, hoek_data$cellFractions.pred)
})
