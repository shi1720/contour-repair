# Literature and contribution record

Search date: 17 October 2026. This is a focused research audit, not a systematic-review claim or a patent search. Absence from these searches does not establish that a result has never appeared elsewhere.

## Screening the author's existing work

The supplied Downloads/research corpus contained 88 PDF files, including venue variants and working dossiers. Text extraction and title/abstract screening searched for contour quadrature, Laurent/complex analysis, analytic continuation, conformal mapping, rational approximation and related terminology. The closest ICERM-labelled paper, **Frobenius Lens**, concerns finite-field cellular-automaton identification using annihilating polynomials. Its problem, mathematics, experiments and code are distinct from the circular analytic recovery problem developed here. Medical contour-volume work concerns planar reconstructions and slice-origin dependence, also a different problem. This screen does not certify the contents of material outside the supplied corpus or unfinished work in other chats.

## Searches and resulting scope decisions

Search families included: contour quadrature with missing/failed evaluations; quadrature erasures; optimal unit-circle and annular quadrature; Poisson kernels and quadrature; incomplete-lattice FFT inference; quadrature updates/downdates; greedy kernel quadrature and submodularity; contour eigensolvers with inexact solves; and one-missing-node exponential convergence. Searches used scholarly web indexing and primary author, journal, arXiv and software sources. Irrelevant hits were not used as evidence.

1. A generic kernel-quadrature proposal was rejected as a novelty claim: optimal circle rules predate this project by decades.
2. Generic greedy retry selection was also rejected as a novel principle: sequential Bayesian quadrature and greedy RKHS functional recovery already exist.
3. The factor-of-two loss from perturbing uniform grids is explicitly documented by Trefethen and Weideman. The paper cites that connection prominently.
4. The chosen contribution is the exact one-erasure minimax expression and constant in the specified geometric Laurent energy space, an independent vanishing-polynomial lower bound, and an erasure-sized FFT recovery and retry specialization.
5. The inverse-block formula is standard algebra and is not claimed as an original matrix identity. General non-submodularity is not new; the four-node exact example is an explicit specialization.
6. No exact duplicate of the stated one-erasure expression was located. Priority remains provisional, especially because the complete text of the 1973 Chawla–Kaul paper was not accessible in this session; its metadata and historical relevance were verified. The manuscript avoids an absolute 'first ever' assertion.

## Primary-source map

| Source | Evidence consulted | Consequence for this paper |
|---|---|---|
| Chawla & Kaul (1973), [DOI](https://doi.org/10.1007/BF01933486) | Bibliographic metadata and historical references; full text unavailable | Optimal unit-circle quadrature is established; no claim to invent it |
| Kress (1978), [DOI](https://doi.org/10.1137/0715027) | Journal abstract; bibliographic pages checked | Error-norm optimality of trapezoidal weights has precedent |
| Trefethen & Weideman (2014), [author PDF](https://people.maths.ox.ac.uk/trefethen/publication/PDF/2014_149.pdf) | Full review, aliasing and uneven-node discussion | Explicitly distinguish the one-erasure exact law from known factor-of-two losses |
| Huszár & Duvenaud (2012), [proceedings PDF](https://auai.org/uai2012/papers/213.pdf) | Full paper, optimal weights, MMD and submodularity discussion | Attribute minimax kernel quadrature and general non-submodularity |
| Karvonen & Särkkä (2018), [author PDF](https://tskarvone.github.io/pdf/KarvonenSarkka2018-fss.pdf) | Full article and journal metadata | Compare complementary inverse-block structure with symmetry reductions |
| Santin, Karvonen & Haasdonk (2022), [preprint](https://arxiv.org/abs/2004.00556) | Abstract and full preprint; final journal metadata | Attribute greedy RKHS functional approximation |
| Stroud, Stein & Lysen (2017), [author PDF](https://jrs390.georgetown.domains/papers/Bayes-MLE-lattice.pdf) | Author version and preprint abstract | Acknowledge incomplete-lattice/circulant computation precedent |
| Gonnet (2010), [preprint](https://arxiv.org/abs/1006.3962) | Preprint abstract and interpolation-downdate description | Removing unavailable function values has prior numerical methods |
| Hale, Higham & Trefethen (2008), [author PDF](https://appliedmaths.sun.ac.za/~nhale/publications/HaleHighamTrefethen2008.pdf) | Full preprint and final title/DOI | Motivation through matrix functions and parallel shifted solves |
| Sakurai & Sugiura (2003), [DOI](https://doi.org/10.1016/S0377-0427(03)00565-X) | Original metadata and official SLEPc reference documentation | Attribute contour eigensolver lineage |
| Imakura, Morikuni & Takayasu (2020), [preprint](https://arxiv.org/abs/1904.06277) | Author abstract and open journal abstract | Distinguish exact mathematical bounds from verified floating-point enclosures |

## Interpretation of results

The strongest finding is a proved information limit within a precisely defined space. The code supplies independent numerical implementations and exact rational checks; experiments cannot prove unqualified novelty. The rational benchmarks show that minimax weights can have worse nominal accuracy than smoother-model Fourier baselines. The manuscript retains those outcomes and the greedy retry gap. No claim of industry-wide superiority or a high numerical acceptance probability is supported.
