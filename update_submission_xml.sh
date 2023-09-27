
#!/bin/bash

# Getting submission file
submission_file=$1

# Path to a temp config file
tmp_sub_xml=/tmp/tmp_submission.xml

# Remove temp submission file if it exists
if [[ -f $tmp_sub_xml ]]
then
	rm $tmp_sub_xml
fi

# Append NEW submission directory in config file 
cat $submission_file | sed -e 's,<Action[0-9]*>,<Action>,g' |  sed -e 's,</Action[0-9]*>,</Action>,g' > $tmp_sub_xml

# Update config file
mv $tmp_sub_xml $submission_file
