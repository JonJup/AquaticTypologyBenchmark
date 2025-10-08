################################################################################
# Script Name:        container_python_test.R
# Description:        Test script. How can I define the path to python in the container. 
#
# Author:             Jonathan Jupke
# Date Created:       2025-09-26
# Last Modified:      2025-09-26
#
# R Version:          R x.x.x
# Required Packages:  package1, package2
#
# Notes:              Any notes or assumptions
################################################################################

python <- Sys.which("python")
print(
        system2(python, "--version")
)