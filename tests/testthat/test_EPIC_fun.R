context("Testing EPIC function")

test_that("Test of bad inputs", {
          a <- matrix(1:20, nrow=5)
          expect_error(EPIC(bulk=a), "There are only 0 signature")
          refStr <- "unknownRef"
          expect_error(EPIC(bulk=melanoma_data$counts, reference=refStr),
                       paste0(refStr, ".* not part of the allowed references"))
          tRef <- BRef; tRef$sigGenes <- NULL
          expect_error(EPIC(bulk=melanoma_data$counts, reference=tRef),
                       "needs to contain .* 'sigGenes'")
          tRef <- BRef;
          tRef$refProfiles.var <- tRef$refProfiles.var[
            nrow(tRef$refProfiles.var):1,]
          expect_error(EPIC(bulk=melanoma_data$counts, reference=tRef),
                       "dimnames of .*refProfiles.*refProfiles.var.*same")
})

test_that("Test for correct result on melanoma data with default input", {
  testFract <- EPIC(melanoma_data$counts)$cellFractions
  expect_equal(testFract, melanoma_data$cellFractions.pred)
})
