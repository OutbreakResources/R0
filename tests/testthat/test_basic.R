context('basic use')

test_that('basic use', {
  threshold<-5e3
  yobs=rep(5,1)#c(10,10)
  pi<-0.5
  res<-R0(yobs,threshold,pi)
  expect_equal(length(res),3)
})