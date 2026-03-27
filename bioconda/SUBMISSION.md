# Bioconda Submission Checklist

Current upstream code release:

- version: `2.0.3`
- source tarball: `https://github.com/kentnf/iTAK/archive/refs/tags/v2.0.3.tar.gz`
- source sha256: `0b6eb8653c3d5e42123d0b671633b3b3b128e8b478621320cc6b644d24fe8308`

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

`Add iTAK 2.0.3`

## Suggested PR Body

```text
This PR adds iTAK 2.0.3 to Bioconda.

Highlights:
- Python package with console entry points `itak` and `iTAK`
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

## Open Item To Confirm Manually

- The repository metadata says `MIT`, but the repository currently does not include a top-level `LICENSE` file. If Bioconda reviewers require `license_file`, add the canonical license text in a separate upstream change first.
