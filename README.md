# cisid-phip-flow-workflow

This workflow is specific to the CISID, using phippery and phip-flow as the base and applied to work as a Terra workflow.

For specific parameters and instructions about phippery and phipflow refer to the original documentation from the Matsen Group.

phippery: https://github.com/matsengrp/phippery
phip-flow: https://github.com/matsengrp/phip-flow

The intention of this repository is to provide a ready-made Dockerfile that contains edgeR, BEER, and Nextflow in order to host and run the phipflow pipeline on Terra.

There were small parameter changes within the main.nf, provided with the Dockerfile & WDL.

This repository is connected to the following dockstore repository: https://dockstore.org/workflows/github.com/klarman-cell-observatory/cisid-phippery-workflow/cisid_phip_flow:main?tab=info to be pulled from/exported to Terra via their builtin interface.

Additional in-house figure generation and downstream analysis to be added in the future on top of base phip-flow.