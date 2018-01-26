context("utility functions")

test_that("binary tree works", {
	b1 <- perfect_binary(1L)
	b3 <- perfect_binary(3L)
	b5 <- perfect_binary(5L)
})

b1 <- perfect_binary(1L)
b3 <- perfect_binary(3L)
b5 <- perfect_binary(5L)

inp <- data.frame(node=1L, rate=1)

n1 <- popsnetwork(b1, inp)
n3 <- popsnetwork(b3, inp)
n5 <- popsnetwork(b5, inp)

test_that("trees are correct and distances are calculated correctly", {
	m1 <- distances_topology(n1)
	m3 <- distances_topology(n3)
	m5 <- distances_topology(n5)

	expect_true(isSymmetric(m1))
	expect_true(isSymmetric(m3))
	expect_true(isSymmetric(m5))

	expect_equal(m1[1, 1], 0)
	expect_equal(m3[1, 1], 0)
	expect_equal(m5[1, 1], 0)

	expect_equal(m1[nrow(m1), nrow(m1)-1], 2)
	expect_equal(m3[nrow(m3), nrow(m3)-7], 6)
	expect_equal(m5[nrow(m5), nrow(m5)-31], 10)

	m1 <- distances_topology(n1, FALSE)
	m3 <- distances_topology(n3, FALSE)
	m5 <- distances_topology(n5, FALSE)

	expect_false(any(m1==-1))
	expect_false(any(m3==-1))
	expect_false(any(m5==-1))

	expect_true(isSymmetric(m1))
	expect_true(isSymmetric(m3))
	expect_true(isSymmetric(m5))
})
