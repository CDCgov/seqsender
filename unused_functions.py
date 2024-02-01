# Check description section for submission.xml
def check_submission_description(config_dict, database):
	# Check if config has description information
	if "Description" not in config_dict.keys():
		print("Error: there is no Submission > " + database + " > Description information in config file.", file=sys.stderr)
		sys.exit(1)
	# Check if description has a title (required)
	try:
		title = config_dict["Description"]["Title"]
	except:
		print("Error: there is no Submission > " + database + " > Description > Title in the config file.", file=sys.stderr)
		sys.exit(1)
	else:
		# Make sure title is not empty or none
		if (title is None) or (title == ""):
			print("Error: Submission > " + database + " > Description > Title in the config file cannot be empty (required).", file=sys.stderr)
			sys.exit(1)
	# Check if description has organization information (required) for NCBI
	if database == "NCBI":
		try:
			organization = config_dict["Description"]["Organization"]
		except:
			print("Error: there is no Submission > " + database + " > Description > Organization in the config file (required).", file=sys.stderr)
			sys.exit(1)
		else:
			# Check organization role
			if (organization["@role"] is None) or (organization["@role"] == ""):
				print("Error: Submission > " + database + " > Description > Organization > @role in the config file cannot be empty.", file=sys.stderr)
				sys.exit(1)
			# Check organization type
			if (organization["@type"] is None) or (organization["@type"] == ""):
				print("Error: Submission > " + database + " > Description > Organization > @type in the config file cannot be empty.", file=sys.stderr)
				sys.exit(1)
			# Check organization name
			if (organization["Name"] is None) or (organization["Name"] == ""):
				print("Error: Submission > " + database + " > Description > Organization > Name in the config file cannot be empty.", file=sys.stderr)
				sys.exit(1)
			# Check contact info
			if (organization["Contact"] is None) or (organization["Contact"] == ""):
				print("Error: Submission > " + database + " > Description > Organization > Contact in the config file cannot be empty.", file=sys.stderr)
				sys.exit(1)
			else:
				contact = organization["Contact"]
				if(contact["@email"] is None) or (contact["@email"] == ""):
					print("Error: Submission > " + database + " > Description > Organization > Contact > @email in the config file cannot be empty.", file=sys.stderr)
					sys.exit(1)
				if(contact["Name"] is None) or (contact["Name"] == ""):
					print("Error: Submission > " + database + " > Description > Organization > Contact > Name in the config file cannot be empty.", file=sys.stderr)
					sys.exit(1)
				else:
					contact_name = contact["Name"]
					if(contact_name["First"] is None) or (contact_name["First"] == ""):
						print("Error: Submission > " + database + " > Description > Organization > Contact > Name > First in the config file cannot be empty.", file=sys.stderr)
						sys.exit(1)
					if(contact_name["Last"] is None) or (contact_name["Last"] == ""):
						print("Error: Submission > " + database + " > Description > Organization > Contact > Name > Last in the config file cannot be empty.", file=sys.stderr)
						sys.exit(1)
	return config_dict["Description"]


# Get authentication token
def get_token(token_file, organism):
	if os.path.exists(token_file) == False:
		print("Error: Authentication token does not exist at: "+token_file, file=sys.stderr)
		print("Error: Either a submission has not been made or token has been moved.", file=sys.stderr)
		print("Please re-authenticate the database to obtain the token again", file=sys.stderr)
		sys.exit(1)
	else:
		with open(token_file) as f:
			token_dict = json.load(f)
			f.close()
		# Extract user credentials
		if type(token_dict) is dict:
			try:
				token_dict = {k.title():v for k,v in token_dict[organism.lower()].items()}
			except:
				print("Error: there is no " + organism.lower() + " information in token file.", file=sys.stderr)
				sys.exit(1)
			else:
				# Check username
				if "Username" not in token_dict.keys():
					print("Error: there is no " + organism.lower() + " > username information in token file.", file=sys.stderr)
					sys.exit(1)
				elif ("Username" in token_dict.keys()) and ((token_dict["Username"] is None) or (token_dict["Username"] == "")):
					print("Error: " + organism + " > username in the token file cannot be empty.", file=sys.stderr)
					sys.exit(1)
				# Check password
				if "Password" not in token_dict.keys():
					print("Error: there is no" + organism.lower() + " > password information in token file.", file=sys.stderr)
					sys.exit(1)
				elif ("Password" in token_dict.keys()) and ((token_dict["Password"] is None) or (token_dict["Password"] == "")):
					print("Error: " + organism.lower() + " > password in the token file cannot be empty.", file=sys.stderr)
					sys.exit(1)
				return token_dict
		else:
			print("Error: Token file is incorrect. Token file must be in a valid json format.", file=sys.stderr)
			sys.exit(1)
			

# Authenticate user credentials
def authenticate(database, organism, config_dict):
	process.check_credentials(config_dict=config_dict, database=database)
	if "GISAID" in database:
		log_file = os.path.join(WORK_DIR, "gisaid_authenticate_logfile.txt")
		# If log file exists, removes it
		if os.path.isfile(log_file) == True:
			os.remove(log_file)
		if "FLU" in organism:
			command = ["python", os.path.join(PROG_DIR, "epiflu_cli/__main__.py"), "authenticate",
				"--token", "gisaid.flu.authtoken", "--user", config_dict["Username"], "--pass", config_dict["Password"],
				"--client_id", config_dict["Client-Id"], "--log", log_file, "--force"]
		elif "COV" in organism:
			command = ["python", os.path.join(PROG_DIR, "epicov_cli/__main__.py"), "authenticate",
				"--token", "gisaid.cov.authtoken", "--user", config_dict["Username"], "--pass", config_dict["Password"],
				"--client_id", config_dict["Client-Id"], "--log", log_file, "--force"]
		# Check status of the command
		proc = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd = WORK_DIR)
		if proc.returncode == 0:
			print("\n" + "User authenticated", file=sys.stdout)
		else:
			print("\n" + "Authentication failed", file=sys.stderr)
			print(proc.stdout, file=sys.stdout)
			print(proc.stderr, file=sys.stdout)
			sys.exit(1)
	elif "NCBI" in database:
		# Login into NCBI FTP Server
		try:
			ftp = ftplib.FTP(NCBI_FTP_HOST)
			ftp.login(user=config_dict["Username"], passwd = config_dict["Password"])
			# Output authenticated message
			print("\n" + "User authenticated", file=sys.stdout)
		except ftplib.all_errors as e:
			print("\n" + "Authentication failed", file=sys.stderr)
			print("\n" + "Error: " + str(e), file=sys.stderr)
			sys.exit(1)

