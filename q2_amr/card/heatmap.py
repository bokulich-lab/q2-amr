import os
import shutil
import subprocess
import tempfile
from distutils.dir_util import copy_tree

import pkg_resources
import q2templates

from q2_amr.card.utils import run_command
from q2_amr.types import CARDAnnotationDirectoryFormat


def heatmap(
    output_dir: str,
    amr_annotation: CARDAnnotationDirectoryFormat,
    clus: str = None,
    cat: str = None,
    display: str = "plain",
    frequency: bool = False,
):
    TEMPLATES = pkg_resources.resource_filename("q2_amr", "assets")
    annotation_dir = str(amr_annotation)
    with tempfile.TemporaryDirectory() as tmp:
        results_dir = os.path.join(tmp, "results")
        json_files_dir = os.path.join(tmp, "json_files")
        os.makedirs(results_dir)
        os.makedirs(json_files_dir)
        for sample in os.listdir(annotation_dir):
            for bin in os.listdir(os.path.join(annotation_dir, sample)):
                for file in os.listdir(os.path.join(annotation_dir, sample, bin)):
                    if file.endswith(".json"):
                        shutil.copy(
                            os.path.join(annotation_dir, sample, bin, file),
                            os.path.join(json_files_dir, f"{sample}_{bin}.json"),
                        )

        run_rgi_heatmap(tmp, json_files_dir, clus, cat, display, frequency)
        change_names(results_dir)
        copy_tree(os.path.join(TEMPLATES, "rgi", "heatmap"), output_dir)
        copy_tree(results_dir, os.path.join(output_dir, "rgi_data"))
    context = {"tabs": [{"title": "Heatmap", "url": "index.html"}]}
    index = os.path.join(TEMPLATES, "rgi", "heatmap", "index.html")
    templates = [index]
    q2templates.render(templates, output_dir, context=context)


class InvalidParameterCombinationError(Exception):
    def __init__(self, message="Invalid parameter combination"):
        self.message = message
        super().__init__(self.message)


def run_rgi_heatmap(tmp, json_files_dir, clus, cat, display, frequency):
    cmd = [
        "rgi",
        "heatmap",
        "--input",
        json_files_dir,
        "--output",
        f"{tmp}/results/heatmap",
        "--display",
        display,
    ]
    if clus:
        cmd.extend(["--clus", clus])
    if cat:
        cmd.extend(["--cat", cat])
    if frequency:
        cmd.append("--frequency")
    if (clus == "both" or clus == "genes") and cat:
        raise InvalidParameterCombinationError(
            "If the parameter clus is set to genes "
            "or both it is not possible to use the "
            "cat parameter"
        )
    try:
        run_command(cmd, tmp, verbose=True)
    except subprocess.CalledProcessError as e:
        raise Exception(
            "An error was encountered while running rgi, "
            f"(return code {e.returncode}), please inspect "
            "stdout and stderr to learn more."
        )


def change_names(results_dir):
    extensions = [".eps", ".csv", ".png"]
    files = os.listdir(results_dir)
    for filename in files:
        if os.path.splitext(filename)[1] in extensions:
            file_ext = os.path.splitext(filename)[1]
            new_filename = "heatmap" + file_ext
            old_path = os.path.join(results_dir, filename)
            new_path = os.path.join(results_dir, new_filename)
            os.rename(old_path, new_path)
