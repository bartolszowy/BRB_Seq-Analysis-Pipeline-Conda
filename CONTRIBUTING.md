# Contributing

We welcome contributions to improve the BRB-seq pipeline!

## Reporting Issues

Found a bug or have a suggestion?

1. Check [existing issues](../../issues)
2. Create a new issue with:
   - Clear description
   - Steps to reproduce (for bugs)
   - Error logs and SLURM job IDs
   - Expected vs actual behavior

## Contributing Code

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-improvement`
3. Make your changes
4. Test on spike-in data
5. Commit: `git commit -m "Add feature: description"`
6. Push: `git push origin feature/my-improvement`
7. Open a Pull Request

## Code Style

**Bash:**
- Use `set -euo pipefail`
- Quote variables: `"${VAR}"`
- Add comments for complex logic

**Python:**
- Follow PEP 8
- Add docstrings
- Use type hints

## Questions?

- [Discussions](../../discussions)
- Lab Slack: #bioinformatics
- Email: [Lab contact]

Thank you! 🧬
