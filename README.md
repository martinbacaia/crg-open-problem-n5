# What is the classical Regge limit for n≥5 particles?

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22820411.svg)](https://doi.org/10.5281/zenodo.22820411)

**Status: preprint, not peer-reviewed.** This note reports a negative/exploratory result and an
open problem, not a resolved theorem. See the disclaimer on the first page of the paper. Feedback,
corrections, and criticism are genuinely welcome — please open an issue.

## What this is

The Classical Regge Growth (CRG) conjecture — that any consistent classical $S$-matrix grows no
faster than $s^2$ in the Regge limit — is well-defined and used to constrain four-particle photon
and graviton $S$-matrices (Chowdhury-Gadde-Gopalka-Halder-Janagal-Minwalla, arXiv:1910.14392).
Extending it to `n≥5` particles first requires generalizing what "the Regge limit" even means, and
this note shows that generalization is genuinely ambiguous — no definition exists in the
literature for `n≥5`, and the closest established framework with a similar name (multi-Regge
kinematics, from gauge-theory multi-particle production) is motivated by different physics.

Using the five-photon operator basis constructed in a companion note (see below), we build three
qualitatively distinct, explicitly validated families of `n=5` Regge trajectories and search for
degree-9 and degree-11 operators compatible with the conjectured bound. The result: candidates
that survive one or two families are killed by the next, and an exact computation of the combined
kernel across all three families, on the complete 16-dimensional degree-11 space, gives exactly
`{0}`. We read this as evidence that the obstruction is in the *definition* of the limit, not in
an incomplete search, and state it as an open problem rather than a resolved negative result.

**Companion note** (the underlying operator basis this note builds on): M. Baca, *The multilinear
gauge-invariant operator basis for five massless photons*, preprint (2026),
https://doi.org/10.5281/zenodo.22819347 — also on GitHub at
https://github.com/martinbacaia/five-photon-operator-basis.

## Structure

```
paper/    — the LaTeX source and compiled PDF (start here)
scripts/  — the Python scripts that produced every number cited in the paper
notes/    — session notes documenting how each result was derived, including bugs found and
            fixed along the way (referenced from the paper's verification-methodology section)
```

## Reproducing the results

Every script is runnable with Python 3 + `sympy` + `mpmath` + `numpy`/`scipy`. There is no single
"run everything" entry point — each script corresponds to a specific claim in the paper. The two
scripts that reproduce the paper's central (negative) result are `grado11_S16_nullspace_combinado.py`
(the dimension-2-to-3 jump from the new atom) and `grado11_S16_nullspace_3familias.py` (the final
combined-kernel computation across all three families, giving rank 16 / kernel `{0}`).

## License

Code (`scripts/`) is MIT-licensed. Text and the manuscript (`paper/`, `notes/`) are licensed under
CC-BY 4.0. See `LICENSE`.

## Provenance

This repository was extracted from a larger, private working repository where this result was
developed alongside other exploratory work; the history here starts fresh from that point rather
than replaying every intermediate commit.
