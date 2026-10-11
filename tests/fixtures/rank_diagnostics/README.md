# Independent diagnostic reference

These exact CSV arrays and expected metrics were produced before production
comparison by the independent scientific reviewer using R posterior 1.6.0.
Generators, package/session provenance, external CmdStan corroboration, and
finite-sample convention differences are kept in the maintainer's local
development records, which are not public; per-file SHA-256 hashes of the
source CSVs and binary arrays are in binary_manifest.json.

`expected.csv` is the unchanged independent output. Frozen/nonfinite original
chains retain PsyRAT's hard-failure rule, even where a raw numerical reference
returns a finite value. Short and one-chain numerical agreement is not proof
of sufficient diagnostic power. Tail ESS concerns 5/95% quantiles, not 95%
interval endpoints. Cauchy fixtures test numerical rank robustness; no
finite-mean population or meaningful mean-MCSE claim follows for Cauchy data.

The executable oracle uses the binary arrays and posterior_binary_reference.csv.
All input doubles are identical in R readBin and MATLAB fread (little endian,
float64, column-major, dimensions in binary_dimensions.csv). Original decimal
fixtures and expectations remain untouched as evidence. R read.csv parsed some
decimals one ULP differently, changing folded median ties in two fixtures.
No diagnostic tolerance or reference implementation changed.
