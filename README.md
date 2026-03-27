# iTAK v2

iTAK is a comprehensive tool designed for the identification and classification of plant transcription factors (TFs), transcriptional regulators (TRs), and protein kinases (PKs) within plant genome sequences. Originally developed in Perl, iTAK v2 has been completely rewritten in Python, offering enhanced performance, usability, and additional features.

## New Features in v2.0.2

- **Rewritten in Python:** Complete rewrite in Python for improved performance and easier maintenance.
- **Extended Database:** The latest Pfam database.
- **Custom Classification:** Users can customize the identification and classification rules of a certain gene family, not limited to TFs and PKs.
- **Bioconda Package:** iTAK v2 is now available as a Bioconda package, simplifying installation and management.

## Installation

iTAK can be easily installed via Bioconda with the following command:

```bash
conda install -c bioconda itak
```

Ensure you have Conda installed and set up before running the installation command. For more detailed installation instructions, including setting up Conda, please refer to the [Bioconda documentation](https://bioconda.github.io/).

For local development, prefer managing the environment with pixi:

```bash
pixi run validate
```

This project expects external tools such as HMMER to be provided by the environment and does not rely on shipping repository-local binaries.

If you prefer not to use pixi, you can still install Python dependencies directly:

```bash
python -m pip install -r requirements.txt
```

Or install the package in editable mode:

```bash
python -m pip install -e .
```

You also need HMMER available on your `PATH`.

## Usage

After installation, iTAK can be run from the command line. Here's a basic example to get you started:

```bash
iTAK.py <sequence_file>
```

From a source checkout, you can run:

```bash
python iTAK.py <sequence_file>
```

You can also use the package-style entry point from the repository root:

```bash
python -m itak <sequence_file>
```

After `pip install -e .`, you can run:

```bash
itak <sequence_file>
```

With pixi, you can run:

```bash
pixi run itak -- <sequence_file>
```

## Testing

Run the validation commands from the repository root:

```bash
bash scripts/validate.sh
```

Or through pixi:

```bash
pixi run validate
```

The smoke test is included in `unittest` discovery and skips automatically when `hmmscan` is not available on `PATH`.

## Documentation

For more information on installation, usage, and customization, please visit the [iTAK documentation page](https://github.com/kentnf/iTAK/wiki).

## Support

If you encounter any issues or have questions regarding iTAK, please open an issue on our [GitHub repository](https://github.com/kentnf/iTAK/issues), and we will be glad to assist.

## Citation

If you use iTAK in your research, please consider citing our publication:

[Zheng Y, Jiao C, Sun H, Rosli HG, Pombo MA, Zhang P, Banf M, Dai X, Martin GB, Giovannoni JJ, Zhao PX, Rhee SY, Fei Z, "iTAK: a program for genome-wide prediction and classification of plant transcription factors, transcriptional regulators, and protein kinases," Molecular Plant, 2016.](https://www.sciencedirect.com/science/article/pii/S1674205216302234)
