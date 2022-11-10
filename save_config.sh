
#!/bin/bash
# Wrapper to run public-repository-submission pipeline

config_file=$1
submission_dir=$2

# Path to a temp config file
tmp_config=/tmp/tmp_config.yaml

# Remove temp config file if it exists
# Because we do not want to keep append to this file
if [[ -f $tmp_config ]]
then
	rm $tmp_config
fi

# Append NEW submission directory in config file 
cat $config_file | awk "NR<3" | grep -v "submission_directory:" >> $tmp_config
echo -e "  submission_directory: $submission_dir" >> $tmp_config
cat $config_file | awk "NR>2" | grep -v "submission_directory:" >> $tmp_config

# Update config file
mv $tmp_config $config_file
