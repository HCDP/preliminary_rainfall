import json, sys, os
from datetime import datetime

config_file = sys.argv[1]
date_s = sys.argv[2]

config = None
# Load the config file specified as the first command line argument
with open(config_file, "r") as f:
    config = json.load(f)

# Inject public access permissions

# Default permission setup for uploading data for public access.
file_permissions = [
    {
        "permission": "READ",
        "recursive": True,
        "username": "public"
    }
]
# inject file_permissions
for upload in config["upload"]:
    upload["file_permissions"] = file_permissions

# Inject auth variables from env
config["agave_options"]["username"] = os.environ["IW_USERNAME"]
config["agave_options"]["password"] = os.environ["IW_PASSWORD"]
config["agave_options"]["api_key"] = os.environ["IW_API_KEY"]
config["agave_options"]["api_secret"] = os.environ["IW_API_SECRET"]

inject_date = datetime.fromisoformat(date_s)
config_s = json.dumps(config, indent = 4)
config_s = data.replace("%y", inject_date.strftime("%Y"))
config_s = data.replace("%m", inject_date.strftime("%m"))
config_s = data.replace("%d", inject_date.strftime("%d"))

# write updated config
with open(config_file, "w") as f:
    f.write(config_s)