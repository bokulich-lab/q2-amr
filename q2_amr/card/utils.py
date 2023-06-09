import json
import subprocess

EXTERNAL_CMD_WARNING = (
    "Running external command line application(s). "
    "This may print messages to stdout and/or stderr.\n"
    "The command(s) being run are below. These commands "
    "cannot be manually re-run as they will depend on "
    "temporary files that no longer exist."
)


def run_command(cmd, cwd, verbose=True):
    if verbose:
        print(EXTERNAL_CMD_WARNING)
        print("\nCommand:", end=" ")
        print(" ".join(cmd), end="\n\n")
    subprocess.run(cmd, check=True, cwd=cwd)


def load_preprocess_card_db(tmp, card_db, operation):
    if operation == "load":
        cmd = ["rgi", "load", "--card_json", str(card_db), "--local"]
    elif operation == "preprocess":
        cmd = ["rgi", "card_annotation", "-i", str(card_db)]
    elif operation == "load_fasta":
        with open(str(card_db)) as f:
            card_data = json.load(f)
            version = card_data["_version"]
        cmd = [
            "rgi",
            "load",
            "-i",
            str(card_db),
            "--card_annotation",
            f"card_database_v{version}.fasta",
            "--local",
        ]
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            f"An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )
