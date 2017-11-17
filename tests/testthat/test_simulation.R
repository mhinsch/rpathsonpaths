context("simulator functions")


el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(1.5, 1, 3))
ext <- data.frame(node=c("A", "B"), rate=c(0.3, 0.1))

test_that("network gets constructed", {
	net <- popsnetwork(el, ext)
	print(net)
	
	expect_equal(edge_list(net)[1:3], el)
	expect_equal(node_list(net)[[2]], c(0.3, 0.1, 0.55, 0.66))
})

net <- popsnetwork(el, ext)
freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)

test_that("setting allele frequencies works", {
	set_allele_freqs(net, list(as.factor(c("A", "C")), freqs))
})

net2 <- set_allele_freqs(net, list(as.factor(c("A", "C")), freqs))

test_that("Dirichlet simulation works", {
	res1 <- popgen_dirichlet(net2, 0.3)
	res2 <- popgen_dirichlet(net2, 0.3, list(as.factor(c("A", "C")), freqs))

	# shouldn't change
	expect_equal(node_list(res1), node_list(res2))
	expect_equal(node_list(net), node_list(res1))
})


test_that("IBM simulation works", {
	# create network
	el <- data.frame(from=c("A", "B", "C"), to=c("C", "C", "D"), rates=c(150, 100, 200))
	ext <- data.frame(node=c("A", "B"), rate=c(300, 100), input=c(1000, 1000))
	net_i <- popsnetwork(el, ext, spread_model="units")

	# set allele frequencies (2 nodes, 3 alleles)
	freqs <- matrix(c(0.1, 0.5, 0.4, 0.9, 0.1, 0), nrow=2, ncol=3, byrow=TRUE)
	ini_freqs <- list(as.factor(c("A", "C")), freqs)

	# or we can initialize and run in one call
	res <- popgen_ibm_mixed(net_i, ini_freqs)
	expect_equal(node_list(res), node_list(net_i))

	# output > input
	expect_error(popgen_ibm_mixed(net, ini_freqs))
	expect_error(popgen_ibm_mixed(net2))
	# no allele frequencies
	expect_error(popgen_ibm_mixed(net_i))
})
