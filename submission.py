from datetime import date

class Submission:
    def __init__(self, cmd_input, file_dict):
        self.submission_name = ""
        self.submission_date = date.today().strftime("%Y-%m-%d")
        self.submission_organism = ""
        self.submission_type = ""
        self.submission_dir = ""
        self.config_dict = ""
        # Database information
        self.
