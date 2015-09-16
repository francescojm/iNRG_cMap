# iNRG_cMap
iterative NetwoRk Guided connectivity Mapping

Drug repositioning, whereby a drug prescribed to treat a given disease is approved to treat a new one, has
gained increased interest in recent years. Various methods attempt to identify repositioning opportunities
by analysing signatures of changes in gene expression induced by drugs, generally as ‘black-box’ tools to
link compounds to transcriptional states in an unbiased manner. We have conceived iNRG_cMap (iterative NetwoRk Guided connectivity Mapping) a strategy that refines these unbiased approaches by making use of prior knowledge about the analysed compounds to compute refined transcriptional signatures of drug response. In this way spurious effects due to non-relevant
secondary drug effects are disentangled and the predictive power of the resulting refined signatures
enhanced. In Iorio et al. (PLoS ONE 2015) we present a proof of principle of this strategy in combination with experimental outcome validations through high-resolution cell biology assays. With this approach we predict that glipizide and
splitomicin perturb microtubules in human cells. In agreement, we find that these two drugs reduce
interphase microtubule growth rates and transiently increase the percentage of mitotic cells. Finally, we
have further verified our findings through the signature reversion paradigm (recently proposed in
computational drug-discovery) on independent dataset derived from a large drug screening study.
This repository contains all the code and data objects needed to reproduce the results presented in our paper and to run iNRG_cMap on user-defined cases.
