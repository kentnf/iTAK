# Bioconda Submission Checklist

Current upstream code release:

- version: `2.0.7`
- source tarball: `https://github.com/kentnf/iTAK/archive/refs/tags/v2.0.7.tar.gz`
- source sha256: `6f02223479a1e07750eeaccd5e7e5c31ed354053d3912edd999a0424b2595ea9`

Current external database release:

- tag: `db-v1`
- database tarball: `https://github.com/kentnf/iTAK/releases/download/db-v1/itak-db-v1.tar.gz`
- database sha256 file: `https://github.com/kentnf/iTAK/releases/download/db-v1/itak-db-v1.tar.gz.sha256`
- database archive sha256: `ca1a2c3057290fddef4769fa1a4dfeef7318a3d1b65eaf17e170495ddca10dd4`

Recipe files in this repository:

- `bioconda/recipe/meta.yaml`
- `bioconda/README.md`

## Submission Steps

1. Fork `bioconda/bioconda-recipes`.
2. Create recipe directory `recipes/itak/`.
3. Copy `bioconda/recipe/meta.yaml` into `recipes/itak/meta.yaml`.
4. Commit the new recipe on your fork.
5. Open a pull request against `bioconda/bioconda-recipes`.

## Suggested PR Title

`Add iTAK 2.0.7`

## Suggested PR Body

```text
This PR adds iTAK 2.0.7 to Bioconda.

Highlights:
- Python package with the `itak` console entry point
- Supports Python 3.8+
- Runtime dependency on `hmmer`
- Large reference database is not bundled in the package
- Database is distributed separately through GitHub Releases and installed by users with `itak db download`

Upstream source:
- https://github.com/kentnf/iTAK

Database release used by the packaged CLI:
- https://github.com/kentnf/iTAK/releases/tag/db-v1

Post-install user workflow:
- `itak db download`
- `itak db verify`
```

## Reviewer Notes

Important packaging detail:

- The repository `database/` directory is intentionally excluded from the package payload.
- The CLI downloads database assets from GitHub Releases on explicit user request.
- This keeps the Bioconda package small and avoids embedding large static data in the recipe source.

## Local Sanity Checks

From this repository, the following already passed before submission:

- `bash scripts/validate.sh`
- `python -m pip wheel . --no-deps`
- `python -m itak db download --path <tmpdir>`
- `python -m itak db verify --path <tmpdir>`

## Packaging Notes Addressed Upstream

- A top-level `LICENSE` file is now included, and the recipe declares `license_file: LICENSE`.
