# iTAK v2

iTAK identifies and classifies plant transcription factors (TFs), transcriptional regulators (TRs), and protein kinases (PKs) from genome or transcriptome sequence data.

## Recommended: pixi

This is the simplest way to run iTAK from this repository because `pixi` manages both Python dependencies and HMMER.

### 1. Install dependencies

```bash
pixi install
```

### 2. Download the database

```bash
pixi run itak -- db download
```

### 3. Run iTAK

```bash
pixi run itak -- <sequence_file>
```

Example:

```bash
pixi run itak -- test_seq
```

## Alternative: local pip install

Use this only if you do not want to use `pixi`.

Requirements:

- Python 3.8+
- `hmmscan` from HMMER available on `PATH`

Install:

```bash
python -m pip install -e .
```

Download the database:

```bash
itak db download
```

Run:

```bash
itak <sequence_file>
```

## Conda / Bioconda

Planned, but not available yet.

The Bioconda recipe has been prepared, but it has not been accepted upstream yet. Do not use the commands below until the Bioconda package is actually published.

Planned install command:

```bash
conda install -c bioconda itak
```

Planned usage:

```bash
itak db download
itak <sequence_file>
```

## Database Location

Inspect the current database location:

```bash
itak db path
```

Inspect the default install target:

```bash
itak db path --target
```

Use a custom database location:

```bash
itak db download --path /path/to/itak-db
ITAK_DB_DIR=/path/to/itak-db itak <sequence_file>
```

Verify database files:

```bash
itak db verify
```

Default lookup/install behavior:

- `ITAK_DB_DIR` if set
- `./database` only when it already contains a valid iTAK database
- `CONDA_PREFIX/share/itak/database` for `pixi`/conda-style environments
- `~/.local/share/itak/database` for user installs

## Validation

```bash
pixi run validate
```

## Documentation

More details are available in the [wiki](https://github.com/kentnf/iTAK/wiki).
