# Contributing to GPlately

Thanks for your interest in GPlately. Bug reports, feature requests and pull requests are all welcome.

## Before you open a pull request

**Discuss first.** For anything beyond a small fix (a typo, a documentation correction, an obvious one-line bug), please open an issue, or comment on an existing one, and get agreement from a maintainer on the approach before opening a pull request.

**Direction.** Pull requests that don't match the agreed direction for the project may be closed without a detailed review. This isn't a judgement on the effort involved; it lets us keep the codebase consistent with limited maintainer time.

For questions about using GPlately, please use [GitHub Discussions](https://github.com/GPlates/gplately/discussions) rather than opening an issue.

## Pull request checklist

- **Link the issue** where the change was discussed.
- **Run the tests** and check they pass: `python -m pytest tests-dir/pytestcases`. Add tests for new functionality.
- **Format the code** with [black](https://black.readthedocs.io/).
- **Save notebooks without outputs.** CI rejects notebooks that contain output cells.
- **Keep existing code working.** Old import paths, function arguments and command-line names that users may rely on should keep working (with a deprecation warning, if they're being replaced).
- **Licensing.** GPlately is released under the [GPL-2.0 license](LICENSE), and contributions are accepted under the same license.
