context("Edgelist utility functions")

edgelist <- data.frame(f=c(0L, 1L, 2L), t=c(2L, 2L, 3L))
edgelistf <- data.frame(f=c("0", "1", "2"), t=c("2", "2", "3"))

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

edgelistc <- data.frame(c(0L, 1L, 2L, 3L), c(2L, 2L, 3L, 0L))
edgelistcf <- data.frame(c("0", "1", "2", "3"), c("2", "2", "3", "0"))
	
test_that("cycles are found", {
	expect_true(cycles(edgelistc))
	expect_true(cycles(edgelistcf))

	c <- cycles(edgelistc, TRUE)
	cf <- cycles(edgelistcf, TRUE)

	expect_equal(sort(c[[1]]), c(0, 2, 3))
	expect_equal(sort(as.character(cf[[1]])), c("0", "2", "3"))
})

edgelistna <- data.frame(c(0L, 1L, 2L, NA), c(2L, 2L, 3L, 3L))
edgelistnaf <- data.frame(c("0", "1", "2", NA), c("2", "2", "3", "3"))

# children

test_that("children are found", {
	expect_equal(children(edgelist), c(list(2L), list(2L), list(3L), list(integer(0))))
})

# parents

test_that("parents are found", {
	expect_equal(parents(edgelist), c(list(integer(0)), list(integer(0)), list(c(0L, 1L)), list(2L)))
})

# descendants

test_that("descendants are found", {
	expect_equal(descendants(edgelist), list(
		data.frame(from=c(0L, 2L), to=c(2L, 3L)),
		data.frame(from=c(1L, 2L), to=c(2L, 3L)),
		data.frame(from=2L, to=3L),
		data.frame()))
})

edgelist_sn <- data.frame(f=c(0L, 1L, 2L, 4L), t=c(2L, 2L, 3L, 5L))

# biggest_subnetwork

test_that("biggest_subnetwork works", {
	expect_equal(biggest_subnetwork(edgelist_sn), edgelist)
})

# depth

test_that("depth works", {
	expect_equal(depth(edgelist, nodes(edgelist)), 2L)
	expect_equal(depth(edgelistf, nodes(edgelist)), 2L)
})

# NA for all

test_that("NA produces errors", {
	expect_error(sources(edgelistna))
	expect_error(sources(edgelistnaf))

	expect_error(sinks(edgelistna))
	expect_error(sinks(edgelistnaf))

	expect_error(cycles(edgelistna))
	expect_error(cycles(edgelistnaf))

	expect_error(depth(edgelistna))
	expect_error(depth(edgelistnaf))

	expect_error(children(edgelistna))
	expect_error(children(edgelistnaf))

	expect_error(parents(edgelistna))
	expect_error(parents(edgelistnaf))

	expect_error(descendants(edgelistna))
	expect_error(descendants(edgelistnaf))
})




