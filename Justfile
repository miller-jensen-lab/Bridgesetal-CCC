# Bridges et al. - Figure reproduction

# Run Fig 6 DBiT-seq spatial analysis (original dataset)
fig-6:
    Rscript Fig6-DBiT-run.R

# Run Fig 6 DBiT-seq spatial analysis (YUMMER 2025 dataset)
fig6-2025:
    Rscript Fig6-DBiT-2025.R

# Install R packages (run once)
install-r:
    Rscript install_packages.R

# Install Python packages (run once)
install-py:
    uv sync
