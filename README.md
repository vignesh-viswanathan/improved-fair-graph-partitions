# Improved Fair Graph Partitions: exact certificate for Theorem 4.3

[![verify certificate](https://github.com/vignesh-viswanathan/improved-fair-graph-partitions/actions/workflows/verify.yml/badge.svg)](https://github.com/vignesh-viswanathan/improved-fair-graph-partitions/actions/workflows/verify.yml)

This repository accompanies the paper

> **Improved Fair Graph Partitions**<br>
> Nisarg Shah (University of Toronto) and Vignesh Viswanathan (University of Massachusetts Amherst)<br>
> [arXiv:2605.03238](https://arxiv.org/abs/2605.03238)

Theorem 4.3 of the paper states that when *k* = 2, every graph with at least 240 nodes admits a balanced partition in the (38/25, 0)-core. Part of the proof is computational. This repository contains that part: an exact certificate (`certificate.json`) and a short, independent checker (`check.py`) that verifies it with exact rational arithmetic. The checker uses no optimizer, numerical tolerance or search, and it does not depend on the program that produced the certificate.

## Running the checker

You need Python 3.9 or later.

```bash
git clone https://github.com/vignesh-viswanathan/improved-fair-graph-partitions.git
cd improved-fair-graph-partitions
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python check.py
```

The check takes about 4–6 minutes on a laptop and should end with

```
PASS: exact certificate for alpha=38/25, augmented N>=240; 157 boxes.
```

The full expected output is in [`verification.txt`](verification.txt). The checks are Python assertions, so `check.py` refuses to run under `python -O`.

## Files

| File | Contents |
| --- | --- |
| [`check.py`](check.py) | The checker (128 lines; needs only SymPy). It rebuilds the constraint rows from the formulas in the paper, then checks every identity, every sign condition and the domain cover. |
| [`certificate.json`](certificate.json) | The certificate: 12 rational identities and a cover of five parameter domains by 157 boxes. It is read only as data. |
| [`requirements.txt`](requirements.txt) | Pinned dependency (`sympy==1.14.0`). |
| [`verification.txt`](verification.txt) | Output of a passing run. |
| [`NOTES.txt`](NOTES.txt) | Technical notes: what the checker proves, the five domains, and the derivation of the A2 sampling coefficient. |

## Scope

The checker replaces the computational step of the proof (Appendix C of the paper). The combinatorial reduction in Appendix B (the two blocking sets, odd-order padding, overlap localization, the exclusion of |B<sub>2</sub>| ≤ 2 and the sampling inequalities) must still be verified by reading the paper. [`NOTES.txt`](NOTES.txt) spells out where the boundary lies.

## Citation

```bibtex
@misc{shah2026improved,
  title         = {Improved Fair Graph Partitions},
  author        = {Shah, Nisarg and Viswanathan, Vignesh},
  year          = {2026},
  eprint        = {2605.03238},
  archivePrefix = {arXiv}
}
```

## License

[MIT](LICENSE)
