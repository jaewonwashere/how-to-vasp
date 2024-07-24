# OCP_screening

Codebase for screening catalysts using the [Open Catalyst Project](https://opencatalystproject.org/index.html). This project is largely based on the Open Catalyst Project tutorials available in Meta AI Research's [fair-chem repository](https://github.com/FAIR-Chem/fairchem).

## Features

- Scripts for determining relevant bulk structures
- Surface cleaving tools
- Adsorption site setting
- Adslab configuration
- Relaxation procedures
- Post-processing utilities

## Setup Guide (UCSB CSC Pod Cluster)

This guide will help you set up your environment on the UCSB cluster for running OCP calculations.

### Prerequisites

- UCSB Center for Scientific Computing account
- Basic understanding of Linux command line interface (CLI)

### Step-by-step Setup

1. **Request an Account**
   - Visit the [UCSB Center for Scientific Computing](https://csc.cnsi.ucsb.edu/) website
   - Allow a few days for account setup

2. **Set Up SSH Client**
   - Choose a client based on your local operating system:
     - Windows: WSL or PuTTY
     - MacOS/Linux: Default terminal
   - Refer to the [UCSB CNSI CSC's HPC overview](https://csc.cnsi.ucsb.edu/sites/default/files/2023-01/HPC_Workshop_Winter_23.pdf) for detailed instructions

3. **Initial Login and Password Change**
   - SSH into your account:
     ```bash
     ssh your-login@pod-login1.cnsi.ucsb.edu
     ```
   - Use the temporary password provided in your email
   - Follow email instructions to change your password

4. **Verify Login**
   - You should see a welcome message:
     ```
     ----------------------------
     
     Welcome to Pod
     For basic documentation to get started please see
     http://csc.cnsi.ucsb.edu/docs/pod-cluster
     ```

5. **Install Miniconda**
   - Run the following commands:
     ```bash
     mkdir -p ~/miniconda3
     wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
     bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
     rm -rf ~/miniconda3/miniconda.sh
     ```
   - Initialize Miniconda:
     ```bash
     ~/miniconda3/bin/conda init bash
     ~/miniconda3/bin/conda init zsh
     ```

6. **Verify Miniconda Installation**
   - Check for the 'miniconda3' folder in your home directory

7. **Copy Environment File**
   - Copy the pre-configured environment file:
     ```bash
     cp /home/jaewon_lee/fairchem_Pod.yml ~
     ```

8. **Create Conda Environment**
   - Run:
     ```bash
     conda env create -f fairchem_Pod.yml
     ```

9. **Activate Environment**
   - Activate the new environment:
     ```bash
     conda activate fair-chem
     ```
   - Your prompt should change to indicate the active environment

10. **Ready to Go!**
    - You're now set up to run OCP calculations

### Useful Linux Commands

| Command | Description | Example |
|---------|-------------|---------|
| `ls` | List files/folders | `ls` |
| `cd` | Change directory | `cd /path/to/directory` |
| `pwd` | Print working directory | `pwd` |

> 💡 **Tip:** The home directory can be accessed with `cd ~` or `cd /home/your_username/`

### Additional Resources

- [CSC's HPC overview](https://csc.cnsi.ucsb.edu/sites/default/files/2023-01/HPC_Workshop_Winter_23.pdf)
- [Linux Command Line Cheat Sheet](https://www.stationx.net/linux-command-line-cheat-sheet/)

For any issues or questions, please contact: jaewon_lee@ucsb.edu
