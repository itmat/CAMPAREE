FROM spliceprep.azurecr.io/spliceprep:1.0.4

RUN apt-get update && apt-get upgrade -y && apt-get install -y git \
# CAMPAREE
 && git clone https://github.com/envisagenics/CAMPAREE.git \
 && git clone https://github.com/itmat/BEERS_UTILS.git \
# The requirements.txt from both repos cause circular dependencies.
# These are the uniqque dependencies from each file
 && poetry add scipy roman pysam plotly matplotlib Sphinx boto3 \
 termcolor prettytable PyYAML /BEERS_UTILS /CAMPAREE

CMD ["bash", "-c", "/docker_entry_point.sh"]