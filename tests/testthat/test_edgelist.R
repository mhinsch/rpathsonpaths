context("Edgelist utility functions")

edgelist <- data.frame(c(0L, 1L, 2L), c(2L, 2L, 3L))
edgelistf <- data.frame(c("0", "1", "2"), c("2", "2", "3"))

test_that("sources are identified correctly", {
	scs <- sources(edgelist)
	scsf <- sources(edgelistf)

	expect_equal(length(scs), 2L)
	expect_equal(length(scsf), 2L)

	expect_equal(scs, c(0, 1))
	expect_equal(scsf, factor(c("0", "1"), levels=c("0", "1", "2", "3")))
})


test_that("sinks are identified correctly", {
	sks <- sinks(edgelist)
	sksf <- sinks(edgelistf)

	expect_equal(length(sks), 1L)
	expect_equal(length(sksf), 1L)

	expect_equal(sks, 3)
	expect_equal(sksf, factor("3", levels=c("0", "1", "2", "3")))
})

	

