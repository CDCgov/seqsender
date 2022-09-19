
#!/bin/bash
# Wrapper to run public-repository-submission pipeline

WORK_DIR=$PWD
RESOURCE_ROOT=/opt/public-repository-submissions
seqsender=$RESOURCE_ROOT/seqsender.py
gisaid_uploader=$RESOURCE_ROOT/gisaid_uploader.py
retrieve=$RESOURCE_ROOT/retrieve.py
check=$RESOURCE_ROOT/check.py
config_file=$RESOURCE_ROOT/config_files/default_config.yaml
metadata_file=$RESOURCE_ROOT/data/metadata.tsv
fasta_file=$RESOURCE_ROOT/data/fasta.fasta

seqsender_usage(){
	echo ""
	echo "Here is a list of commands:"
	echo ""

	echo "To create the files for submission and adds to automated submission pipeline and starts submission process:"
	echo -e "\t Usage: seqsender submit [-h] --unique_name <> --fasta <> --metadata <> [--config <>] [--test] [--overwrite]"
	echo ""
	
	echo "To create the files for submission:"
	echo -e "\t Usage: seqsender prep [-h] --unique_name <> --fasta <> --metadata <> [--config <>] [--test] [--overwrite]"
	echo ""
	
	echo "To update process of all sequences in submission pipeline, performs submission to subsequent databases based on submission status:"
	echo -e "\t Usage: seqsender update_submissions [-h]"
	echo ""

	echo "To output BioProject ID needed to perform test submissions:"
	echo -e "\t Usage: seqsender test_bioproject [-h] --unique_name <> [--config <>]"
	echo ""

	echo "To perform manual submission to GISAID:"
	echo -e "\t Usage: seqsender gisaid [-h] --unique_name <> [--config <>] [--test] [--overwrite]"
	echo ""

	echo "To perform manual submission to Genbank:"
	echo -e "\t Usage: seqsender genbank [-h] --unique_name <> [--config <>] [--test] [--overwrite]"
	echo ""

	echo "To perform manual submission to BioSample:"
	echo -e "\t Usage: seqsender biosample [-h] --unique_name <> [--config <>] [--test] [--overwrite]"
	echo ""

	echo "To perform manual submission to SRA:"
	echo -e "\t Usage: seqsender sra [-h] --unique_name <> [--config <>] [--test] [--overwrite]"
	echo ""

	echo "To perform manual submission to BioSample/SRA:"
	echo -e "\t Usage: seqsender biosample_sra [-h] --unique_name <> [--config <>] [--test] [--overwrite]"
}

gisaid_uploader_usage(){
	echo ""
	echo "To authenticate GISAID submissions:"
	echo -e "\t Usage: gisaid_uploader [-h] [--proxy <>] [--debug] [--version] [-a AUTHFILE] [-t AUTHTOKEN] [-l LOGFILE] ctx {authenticate,upload,revoke} ..."
}

retrieve_usage(){
	echo ""
	echo "To retrieve a sample file used by seqsender to submit submissions such as config file, metadata file, and required columns file."
	echo -e "\t Usage: retrieve --sample_file <> [--outdir <>]"
}

check_usage(){
	echo ""
	echo "To check status of a submission:"
	echo -e "\t Usage: check --unique_name <>"
}

required_flags(){
	echo ""
	echo "Required flags:"
	echo -e "\t --unique_name <> \t Unique identifier"
	echo -e "\t --metadata <> \t\t Metadata file"
	echo -e "\t --fasta <> \t\t Fasta file"
	echo -e "\t --sample_file <> \t Name of a file to retrieve: config/metadata/required_columns"
	echo ""
}

optional_flags(){
	echo ""
	echo "Optional flags:"
	echo -e "\t -h, --help \t Show this help message and exit."
	echo -e "\t --config <> \t If using a different config file than the default config. Provide the full name of the config file stored in config_files folder."
	echo -e "\t --test \t Performs test submission to NCBI. Does not perform test submission to GISAID. You must used authenticated CID for test submission to GISAID."
	echo -e "\t --overwrite \t Overwrites an existing submission on NCBI FTP. Used to update errored submissions."
	echo -e "\t --outdir \t Output directory to return the file"
	echo ""
}

usage(){
	seqsender_usage;
	gisaid_uploader_usage;
	retrieve_usage;
	check_usage;
	required_flags;
	optional_flags;
}

# Check commands that user provided
case $1 in 
	seqsender) api=$seqsender;;
	gisaid_uploader) api=$gisaid_uploader;;
	retrieve) api=$retrieve;;
	check) api=$check;;
	-h | --help | "") usage; exit 1;;
	*) echo -e "\n$1 is not a valid command."; usage; exit 1;;
esac

# Check if upload_log.csv exists. If exists, copy it to pipeline directory
if [[ -f $WORK_DIR/upload_log.csv ]]
then
	cp $WORK_DIR/upload_log.csv $RESOURCE_ROOT/upload_log.csv 2>/dev/null
	
	# Change directory and config paths in upload_log.csv to docker directory
	python $RESOURCE_ROOT/extract_upload_log.py --log $RESOURCE_ROOT/upload_log.csv
	
	# Check python exit code
	if [[ $? != 0 ]]
	then 
		exit 1
	fi
fi

# Check if required_columns.yaml exists. If exists, copy it to pipeline directory
if [[ -f $WORK_DIR/required_columns.yaml ]]
then
    cp $WORK_DIR/required_columns.yaml $RESOURCE_ROOT/config_files/required_columns.yaml 2>/dev/null
fi

# Copy submit.ready to pipeline directory. If exists, copy it to pipeline directory
if [[ -f $WORK_DIR/submit.ready ]]
then
    cp $WORK_DIR/submit.ready $RESOURCE_ROOT/submit.ready 2>/dev/null
fi

# Copy gisaid_uploader.log to pipeline directory. If exists, copy it to pipeline directory
if [[ -f $WORK_DIR/gisaid_uploader.log ]]
then
    cp $WORK_DIR/gisaid_uploader.log $RESOURCE_ROOT/gisaid_uploader.log 2>/dev/null
fi

# Copy gisaid_uploader.authtoken to pipeline directory. If exists, copy it to pipeline directory
if [[ -f $WORK_DIR/gisaid_uploader.authtoken ]]
then
    cp $WORK_DIR/gisaid_uploader.authtoken $RESOURCE_ROOT/gisaid_uploader.authtoken 2>/dev/null
fi

# Check commands that user provided
i=2
n=$#
arg_list=""

while [[ i -le $n ]];
do
	var=${@:$i:1}
	#echo $i
	#echo $var
	case "$var" in
		--config)
			# Get arg position
			config_pos=$i
			# Get arg value
			config_path=`echo "${@:($config_pos+1):1}" | sed "s,^./,,g"`
			# echo $config_path
			# Check if config value is empty
			case $config_path in 
				"") echo "Error: config file not provided"; exit 1;;
				*) 
					# Check if config file exists. If exists, copy it to pipeline directory
					provided_config=`echo "$WORK_DIR/$config_path" | sed "s,//*,/,g"`
					config_name=`echo "$provided_config" | awk -F "/" '{print $NF}'`
					config_file=`echo "$RESOURCE_ROOT/config_files/$config_name" | sed "s,//*,/,g"`
					# echo $provided_config
					# echo $config_file
					if [[ -f $provided_config ]]
					then
						# Change submission directory in config files to docker directory
						python $RESOURCE_ROOT/extract_config.py --config $provided_config
						
						# Check python exit code
						if [[ $? != 0 ]]
						then 
							exit 1
						fi

						arg_list+="--config $config_name "
						i=$(($i+2))
					else
						echo "Error: Config file does not exist at: ./$config_path"
						exit 1
					fi
					;;
			esac
			;;
		--metadata) 
			# Get arg position
			metadata_pos=$i
			# Get arg value
			metadata_path=`echo "${@:($metadata_pos+1):1}" | sed "s,^./,,g"`
			# echo $metadata_path
			# Check if metadata value is empty
			case $metadata_path in 
				"") echo "Error: metadata file not provided"; exit 1;;
				*)	
					# Check if metadata file exists. If exists, copy it to pipeline directory
					provided_metadata=`echo "$WORK_DIR/$metadata_path" | sed "s,//*,/,g"`
					# echo $provided_metadata
					if [[ -f $provided_metadata ]]
					then
						# Check if sra_file_path_1 and sra_file_path_2 provided in metadata file exists
						python $RESOURCE_ROOT/extract_metadata.py --metadata $provided_metadata 
						
						# Check python exit code
						if [[ $? != 0 ]]
						then 
							exit 1
						fi

						arg_list+="--metadata $metadata_file "
						i=$(($i+2))
					else
						echo "Error: metadata file does not exist at: ./$metadata_path"
						exit 1
					fi
					;;
			esac
			;; 
		--fasta) 
			# Get arg position
			fasta_pos=$i
			# Get arg value
			fasta_path=`echo "${@:($fasta_pos+1):1}" | sed "s,^./,,g"`
			# echo $fasta_path
			# Check if fasta value is empty
			case $fasta_path in 
				"") echo "Error: fasta file not provided"; exit 1;;
				*)
					# Check if fasta file exists. If exists, copy it to pipeline directory
					provided_fasta=`echo "$WORK_DIR/$fasta_path" | sed "s,//*,/,g"`
					# echo $provided_fasta
					if [[ -f $provided_fasta ]]
					then
						# Copy file to pipeline directory
						cp $provided_fasta $fasta_file 2>/dev/null
						arg_list+="--fasta $fasta_file "
						i=$(($i+2))
					else
						echo "Error: fasta file does not exist at: ./$fasta_path"
						exit 1
					fi
					;;
			esac
			;; 
		--unique_name)
			# Get arg position
			unique_name_pos=$i
			# Get arg value
			unique_name=${@:($unique_name_pos+1):1}
			arg_list+="--unique_name $unique_name "
			i=$(($i+2))
			;; 
		*)
			arg_pos=$i
			arg_val=${@:$arg_pos:1}
			arg_list+="$arg_val "
			i=$(($i+1))
			;; 	
	esac
done

# Execute the submission pipeline
python $api $arg_list

# Check python exit code
if [[ $? != 0 ]]
then 
	exit 1
fi

# Check if upload_log.csv exists. If exists, copy it back to user directory
if [[ -f $RESOURCE_ROOT/upload_log.csv ]]
then
	# Change directory and config paths in upload_log.csv to user directory
	python $RESOURCE_ROOT/update_upload_log.py --log $RESOURCE_ROOT/upload_log.csv
	
	# Check python exit code
	if [[ $? != 0 ]]
	then 
		exit 1
	fi

	cp $RESOURCE_ROOT/upload_log.csv $WORK_DIR/upload_log.csv 2>/dev/null
fi

# Change directory in config file to user directory
case $1 in 
	seqsender)
		case $2 in	
			submit | prep | test_bioproject | gisaid | genbank | biosample | sra | biosample_sra)

				case $3 in
					"" | -h | --help) ;;
					*)
						# Change directory in config file to user directory
						python $RESOURCE_ROOT/update_config.py --config $config_file --unique_name $unique_name
						
						# Check python exit code
						if [[ $? != 0 ]]
						then 
							exit 1
						fi
						;;
				esac
				;;

			*) ;;
		esac
		;;
	*) ;;
esac

# Copy submit.ready back to working directory so users can access it
if [[ -f $RESOURCE_ROOT/submit.ready ]]
then
    cp $RESOURCE_ROOT/submit.ready $WORK_DIR/submit.ready 2>/dev/null
fi

# Copy gisaid_uploader.log back to working directory so users can access it
if [[ -f $RESOURCE_ROOT/gisaid_uploader.log ]]
then
    cp $RESOURCE_ROOT/gisaid_uploader.log $WORK_DIR/gisaid_uploader.log 2>/dev/null
fi

# Copy gisaid_uploader.authtoken back to working directory so users can access it
if [[ -f $RESOURCE_ROOT/gisaid_uploader.authtoken ]]
then
    cp $RESOURCE_ROOT/gisaid_uploader.authtoken $WORK_DIR/gisaid_uploader.authtoken 2>/dev/null
fi
