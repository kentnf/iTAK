# iTAK v2

iTAK identifies and classifies plant transcription factors (TFs), transcriptional regulators (TRs), and protein kinases (PKs) from genome or transcriptome sequence data. The current codebase is the Python rewrite of the original Perl implementation.

## Installation

For local development, prefer `pixi` so Python dependencies and HMMER are managed together:

```bash
pixi install
pixi run validate
```

If you do not want to use `pixi`, install the package and runtime dependencies yourself:

```bash
python -m pip install -e .
```

iTAK expects `hmmscan` from HMMER to be available on `PATH`.

## Database Management

The iTAK database is no longer bundled into the Python package. This keeps the package small enough for package managers such as Bioconda and makes database updates independent from code releases.

Inspect the current database location:

```bash
python -m itak db path
```

Inspect the default install target:

```bash
python -m itak db path --target
```

Download and install the database from GitHub Releases:

```bash
python -m itak db download
```

To build the database release asset from this repository checkout:

```bash
bash scripts/build_db_release.sh db-v1
```

Install into a custom location:

```bash
python -m itak db download --path /path/to/itak-db
```

Verify the required database files:

```bash
python -m itak db verify
```

To force a custom database location at runtime, set `ITAK_DB_DIR`.

## Usage

Run the classic script entry point:

```bash
python iTAK.py <sequence_file>
```

Run the package entry point:

```bash
python -m itak <sequence_file>
```

After installation, the console scripts are also available:

```bash
itak <sequence_file>
iTAK <sequence_file>
```

With `pixi`:

```bash
pixi run itak -- <sequence_file>
pixi run itak -- db download
```

## Testing

Run the validation suite from the repository root:

```bash
bash scripts/validate.sh
```

Or through `pixi`:

```bash
pixi run validate
```

The smoke test is included in `unittest` discovery and skips automatically when `hmmscan` is not available on `PATH`.

## Bioconda Packaging

The recommended distribution split is:

- Bioconda package: Python code plus runtime dependency on HMMER.
- GitHub Releases: database archive and optional checksum files.
- User post-install step: `itak db download`.

A starter Bioconda recipe template is included under `bioconda/recipe/`.

## Documentation

For more information on installation, usage, and customization, see the [iTAK wiki](https://github.com/kentnf/iTAK/wiki).

## Support

Open an issue on the [GitHub repository](https://github.com/kentnf/iTAK/issues) if you hit a bug or packaging problem.

## Citation

If you use iTAK in your research, please cite:

[Zheng Y, Jiao C, Sun H, Rosli HG, Pombo MA, Zhang P, Banf M, Dai X, Martin GB, Giovannoni JJ, Zhao PX, Rhee SY, Fei Z, "iTAK: a program for genome-wide prediction and classification of plant transcription factors, transcriptional regulators, and protein kinases," Molecular Plant, 2016.](https://www.sciencedirect.com/science/article/pii/S1674205216302234)
