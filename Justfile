# Bridges et al. - Figure reproduction

# Run Fig 6 DBiT-seq spatial analysis
fig-6:
    Rscript Fig6-DBiT-run.R

# Install R packages (works with Homebrew R)
install-r:
    Rscript install_packages.R

# Install R packages with rv (requires R Framework, not Homebrew R)
# See: https://github.com/a2-ai/rv/issues/XXX
install-r-rv:
    rv sync

# Install Python packages (run once)
install-py:
    uv sync
