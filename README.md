# iTAK v2

iTAK is a comprehensive tool designed for the identification and classification of plant transcription factors (TFs), transcriptional regulators (TRs), and protein kinases (PKs) within plant genome sequences. Originally developed in Perl, iTAK v2 has been completely rewritten in Python, offering enhanced performance, usability, and additional features.

## New Features in v2.0.1

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

## Usage

After installation, iTAK can be run from the command line. Here's a basic example to get you started:

```bash
iTAK.py <sequence_file>
```

## Documentation

For more information on installation, usage, and customization, please visit the [iTAK documentation page](https://github.com/kentnf/iTAK/wiki).

## Support

If you encounter any issues or have questions regarding iTAK, please open an issue on our [GitHub repository](https://github.com/kentnf/iTAK/issues), and we will be glad to assist.

## Citation

If you use iTAK in your research, please consider citing our publication:

[Zheng Y, Jiao C, Sun H, Rosli HG, Pombo MA, Zhang P, Banf M, Dai X, Martin GB, Giovannoni JJ, Zhao PX, Rhee SY, Fei Z, "iTAK: a program for genome-wide prediction and classification of plant transcription factors, transcriptional regulators, and protein kinases," Molecular Plant, 2016.](https://www.sciencedirect.com/science/article/pii/S1674205216302234)
