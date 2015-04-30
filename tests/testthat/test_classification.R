# library(randomForest)
library(parallelRandomForest)
# library(testthat)

aq <- airquality[ which( apply(airquality, 1, function(row) sum(is.na(row)) ) == 0), ]

ozone <- factor(aq$Ozone > median(aq$Ozone))
aq$Solar.R  <- aq$Solar.R / 2
aq$Wind  <- 10*aq$Wind


x_raw <- as.matrix( aq[,2:6] )
storage.mode(x_raw) <- "raw"
set.seed(131)
rf_raw <- randomForest(x=x_raw, y=ozone, xtest=x_raw, mtry=3, importance=F, na.action=na.omit, keep.forest=T)

x <- as.matrix( aq[,2:6] )
set.seed(131)
rf_double <- randomForest(x=x, y=ozone, xtest=x, mtry=3, importance=F, na.action=na.omit, keep.forest=T)

# dput(rf_raw$predicted)
# dput(rf_double$predicted)
# print(rf_raw)

expected_prediction_double <- structure(c(2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                          1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 1L, 1L, 2L, 2L, 1L, 1L,
                                          1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 1L,
                                          2L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L,
                                          2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 1L, 2L, 2L,
                                          2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L,
                                          1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L, 1L, 2L, 1L, 1L, 1L, 1L,
                                          1L, 1L), .Label = c("FALSE", "TRUE"), class = "factor", .Names = c("1",
                                           "2", "3", "4", "7", "8", "9", "12", "13", "14", "15", "16", "17",
                                           "18", "19", "20", "21", "22", "23", "24", "28", "29", "30", "31",
                                           "38", "40", "41", "44", "47", "48", "49", "50", "51", "62", "63",
                                           "64", "66", "67", "68", "69", "70", "71", "73", "74", "76", "77",
                                           "78", "79", "80", "81", "82", "85", "86", "87", "88", "89", "90",
                                           "91", "92", "93", "94", "95", "99", "100", "101", "104", "105",
                                           "106", "108", "109", "110", "111", "112", "113", "114", "116",
                                           "117", "118", "120", "121", "122", "123", "124", "125", "126",
                                           "127", "128", "129", "130", "131", "132", "133", "134", "135",
                                           "136", "137", "138", "139", "140", "141", "142", "143", "144",
                                           "145", "146", "147", "148", "149", "151", "152", "153"))

expected_prediction_raw <- structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 1L,
                                       1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                                       2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                                       2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L,
                                       1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                                       1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L,
                                       2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 2L,
                                       1L, 1L), .Label = c("FALSE", "TRUE"), class = "factor",
                                     .Names = c("1",
                                      "2", "3", "4", "7", "8", "9", "12", "13", "14", "15", "16", "17",
                                      "18", "19", "20", "21", "22", "23", "24", "28", "29", "30", "31",
                                      "38", "40", "41", "44", "47", "48", "49", "50", "51", "62", "63",
                                      "64", "66", "67", "68", "69", "70", "71", "73", "74", "76", "77",
                                      "78", "79", "80", "81", "82", "85", "86", "87", "88", "89", "90",
                                      "91", "92", "93", "94", "95", "99", "100", "101", "104", "105",
                                      "106", "108", "109", "110", "111", "112", "113", "114", "116",
                                      "117", "118", "120", "121", "122", "123", "124", "125", "126",
                                      "127", "128", "129", "130", "131", "132", "133", "134", "135",
                                      "136", "137", "138", "139", "140", "141", "142", "143", "144",
                                      "145", "146", "147", "148", "149", "151", "152", "153"))


test_that("Predicted classification values", {
  # not stable!
  # expect_equal(expected_prediction_raw, rf_raw$predicted);
  expect_equal(expected_prediction_double, rf_double$predicted);
} )


