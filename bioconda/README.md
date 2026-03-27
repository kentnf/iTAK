## Bioconda Packaging Notes

Use the recipe in `bioconda/recipe/meta.yaml` as the submission base for `bioconda-recipes`.

Recommended distribution split:

- Conda package: Python code plus runtime dependency on `hmmer`
- GitHub Releases: database tarball and matching `.sha256`
- User setup step after installation: `itak db download`

This keeps the Bioconda package small and avoids bundling the large `database/` directory into the recipe source.

Suggested release assets for the database:

- `itak-db-v1.tar.gz`
- `itak-db-v1.tar.gz.sha256`

Suggested user workflow after `conda install -c bioconda itak`:

```bash
itak db download
itak db verify
```

Before submitting a new Bioconda update:

1. Update the recipe `version`.
2. Update the source `sha256`.
3. Confirm `itak --help` and `python -m itak db path --target` work in a clean conda environment.
