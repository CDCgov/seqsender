
#!/bin/bash
# Wrapper to run public-repository-submission pipeline

config_file=$1
submission_directory=$2
outfile=$3

# Path to a temp config file
tmp_config=/tmp/tmp_config.yaml

# Remove temp config file if it exists
if [[ -f $tmp_config ]]
then
	rm $tmp_config
fi

# Append the submission directory to config file using docker directory
cat $config_file | awk "NR<3" | grep -v "submission_directory:" >> $tmp_config
echo -e "  submission_directory: $submission_directory" >> $tmp_config
cat $config_file | awk "NR>2" | grep -v "submission_directory:" >> $tmp_config

# Save a copy back to users
mv $tmp_config $outfile