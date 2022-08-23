###############################################################
### THIS IS NOT REQUIRED ANYMORE SINCE INSTALLING VIA CONDA ###
###############################################################

# Setting up ichorCNA using original repo while still keeping our own snakefiles.
# 1. cd into workflow/scripts/2_copy_number_analysis/.
# 2. Clone ichorCNA.
# git clone --depth=1 --branch=master https://github.com/broadinstitute/ichorCNA.git
# 3. Remove the .git folder.
# rm -IvRf ichorCNA/.git/
# 4. Install ichorCNA.
# R CMD INSTALL ichorCNA

# After returning to workflow/rules/2_copy_number_analysis/ you can start using ichorCNA.
