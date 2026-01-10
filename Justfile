# Bridges et al. - Figure reproduction

# Run Fig 6 DBiT-seq spatial analysis
fig-6:
    Rscript Fig6-DBiT-run.R

# Install R packages (run once)
install-r:
    Rscript install_packages.R

# Install Python packages (run once)
install-py:
    uv sync
