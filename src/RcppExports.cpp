// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "rpathsonpaths_types.h"
#include <Rcpp.h>

using namespace Rcpp;

// sources
IntegerVector sources(const DataFrame& edge_list);
RcppExport SEXP rpathsonpaths_sources(SEXP edge_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type edge_list(edge_listSEXP);
    rcpp_result_gen = Rcpp::wrap(sources(edge_list));
    return rcpp_result_gen;
END_RCPP
}
// colour_network
IntegerVector colour_network(const DataFrame& edge_list);
RcppExport SEXP rpathsonpaths_colour_network(SEXP edge_listSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type edge_list(edge_listSEXP);
    rcpp_result_gen = Rcpp::wrap(colour_network(edge_list));
    return rcpp_result_gen;
END_RCPP
}
// popsnetwork
XPtr<Net_t> popsnetwork(const DataFrame& links, const DataFrame& external, double transmission);
RcppExport SEXP rpathsonpaths_popsnetwork(SEXP linksSEXP, SEXP externalSEXP, SEXP transmissionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const DataFrame& >::type links(linksSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type external(externalSEXP);
    Rcpp::traits::input_parameter< double >::type transmission(transmissionSEXP);
    rcpp_result_gen = Rcpp::wrap(popsnetwork(links, external, transmission));
    return rcpp_result_gen;
END_RCPP
}
// print_popsnetwork
void print_popsnetwork(const XPtr<Net_t>& p_net);
RcppExport SEXP rpathsonpaths_print_popsnetwork(SEXP p_netSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const XPtr<Net_t>& >::type p_net(p_netSEXP);
    print_popsnetwork(p_net);
    return R_NilValue;
END_RCPP
}
// print_popsnode
void print_popsnode(const XPtr<Node_t>& p_node);
RcppExport SEXP rpathsonpaths_print_popsnode(SEXP p_nodeSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const XPtr<Node_t>& >::type p_node(p_nodeSEXP);
    print_popsnode(p_node);
    return R_NilValue;
END_RCPP
}
// set_allele_freqs
XPtr<Net_t> set_allele_freqs(const XPtr<Net_t>& p_net, const List& ini_dist);
RcppExport SEXP rpathsonpaths_set_allele_freqs(SEXP p_netSEXP, SEXP ini_distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const XPtr<Net_t>& >::type p_net(p_netSEXP);
    Rcpp::traits::input_parameter< const List& >::type ini_dist(ini_distSEXP);
    rcpp_result_gen = Rcpp::wrap(set_allele_freqs(p_net, ini_dist));
    return rcpp_result_gen;
END_RCPP
}
// spread_dirichlet
XPtr<Net_t> spread_dirichlet(const XPtr<Net_t>& p_net, double theta, Nullable<List> ini_dist);
RcppExport SEXP rpathsonpaths_spread_dirichlet(SEXP p_netSEXP, SEXP thetaSEXP, SEXP ini_distSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const XPtr<Net_t>& >::type p_net(p_netSEXP);
    Rcpp::traits::input_parameter< double >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< Nullable<List> >::type ini_dist(ini_distSEXP);
    rcpp_result_gen = Rcpp::wrap(spread_dirichlet(p_net, theta, ini_dist));
    return rcpp_result_gen;
END_RCPP
}
// get_popsnode
XPtr<Node_t> get_popsnode(const XPtr<Net_t>& p_net, SEXP id);
RcppExport SEXP rpathsonpaths_get_popsnode(SEXP p_netSEXP, SEXP idSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const XPtr<Net_t>& >::type p_net(p_netSEXP);
    Rcpp::traits::input_parameter< SEXP >::type id(idSEXP);
    rcpp_result_gen = Rcpp::wrap(get_popsnode(p_net, id));
    return rcpp_result_gen;
END_RCPP
}
// draw_isolates_popsnode
IntegerVector draw_isolates_popsnode(const XPtr<Node_t>& p_node, int n);
RcppExport SEXP rpathsonpaths_draw_isolates_popsnode(SEXP p_nodeSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const XPtr<Node_t>& >::type p_node(p_nodeSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(draw_isolates_popsnode(p_node, n));
    return rcpp_result_gen;
END_RCPP
}
// draw_isolates_popsnetwork
DataFrame draw_isolates_popsnetwork(const XPtr<Net_t>& p_net, const DataFrame& samples);
RcppExport SEXP rpathsonpaths_draw_isolates_popsnetwork(SEXP p_netSEXP, SEXP samplesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const XPtr<Net_t>& >::type p_net(p_netSEXP);
    Rcpp::traits::input_parameter< const DataFrame& >::type samples(samplesSEXP);
    rcpp_result_gen = Rcpp::wrap(draw_isolates_popsnetwork(p_net, samples));
    return rcpp_result_gen;
END_RCPP
}
// edge_list
DataFrame edge_list(const XPtr<Net_t>& p_net);
RcppExport SEXP rpathsonpaths_edge_list(SEXP p_netSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const XPtr<Net_t>& >::type p_net(p_netSEXP);
    rcpp_result_gen = Rcpp::wrap(edge_list(p_net));
    return rcpp_result_gen;
END_RCPP
}
// node_list
DataFrame node_list(const XPtr<Net_t>& p_net);
RcppExport SEXP rpathsonpaths_node_list(SEXP p_netSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const XPtr<Net_t>& >::type p_net(p_netSEXP);
    rcpp_result_gen = Rcpp::wrap(node_list(p_net));
    return rcpp_result_gen;
END_RCPP
}
