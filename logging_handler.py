
# Python Libraries
from loguru import logger
import sys
# --------------------------------------------------
# Setup logger output
# --------------------------------------------------
logger.remove(0) # Remove default logger outputs REQUIRED
# This allows new outputs for default logger level formats
cmd_line_format = "<level>{level} | {message}</level>"
cmd_line_error_format = "<level>{level} | {module}:{function}:{line} - {message}</level>"
# --------------------------------------------------
# Setup HEADER
# Display SeqSender header output for run info
# --------------------------------------------------
logger.level("HEADER", no=9, color = "<light-cyan>")
logger.add("seqsender.debug.log", format = "<level>{message}</level>", level = "HEADER", colorize = False, filter = lambda record: record["level"].name == "HEADER")
# --------------------------------------------------
# Setup DEBUG
# Record helpful information
# --------------------------------------------------
logger.add("seqsender.debug.log", level = "DEBUG", colorize = False)
# --------------------------------------------------
# Setup INFO
# Record what SeqSender is doing
# --------------------------------------------------
logger.add("seqsender.log", level = "INFO", colorize = False)
# --------------------------------------------------
# Setup SUCCESS
# Record what steps SeqSender has finished
# --------------------------------------------------
# logger.add(sys.stdout, format = cmd_line_format, level = "SUCCESS", colorize = True)
logger.add("seqsender.log", level = "SUCCESS", colorize = False)
# --------------------------------------------------
# Setup WARNING
# Record potential errors SeqSender identified but did not stop
# --------------------------------------------------
# logger.add(sys.stderr, format = cmd_line_format, level = "WARNING", colorize = True)
logger.add("seqsender.error.log", level = "WARNING", colorize = False)
# --------------------------------------------------
# Setup ERROR
# Record error SeqSender caught and identified which caused it to stop
# --------------------------------------------------
# logger.add(sys.stderr, format = cmd_line_error_format, level = "ERROR", colorize = True)
logger.add("seqsender.error.log", level = "ERROR", colorize = False)
# --------------------------------------------------
# Setup CRITICAL
# Record error SeqSender could not identify (i.e. code bugs)
# --------------------------------------------------
# logger.add(sys.stderr, format = cmd_line_error_format, level = "CRITICAL", colorize = True)
logger.add("seqsender.error.log", level = "CRITICAL", colorize = False)

CONFIGURED_LOGGER = logger
