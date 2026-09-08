test_that("XGBoost default threads and train API work through fit_model", {
  skip_if_not_installed("xgboost")
  set.seed(84)
  df <- data.frame(x = rnorm(100), z = rnorm(100))
  df$y <- 2 * df$x + rnorm(100, sd = .1)
  expect_warning(fit <- fit_model(y ~ .,
    data = df, method = "xgb", nrounds = 30,
    params = list(max_depth = 2, eta = .2, objective = "reg:squarederror", seed = 84)
  ), NA)
  expect_gt(fit$metrics$r_squared, .8)
  expect_true(all(is.finite(predict(fit$model, stats::model.matrix(y ~ ., df)))))
  expect_named(fit$coefficients, c("variable", "gain", "cover", "frequency"))
  expect_warning(one <- fit_xgb(y ~ ., df, nrounds = 5, nthread = 1), NA)
  expect_true(is.finite(one$metrics$r_squared))
  expect_error(fit_xgb(y ~ ., df, nthread = -2), "nthread must")
})
