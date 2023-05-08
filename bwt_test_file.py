from qiime2 import Artifact
from qiime2.plugins import amr

reads = Artifact.load("/Users/vinzent/Desktop/bokulich_project/data/amr_bwt/input/reads.qza")
db = Artifact.load("/Users/vinzent/Desktop/bokulich_project/data/amr_bwt/input/card_db.qza")

amr.visualizers.bwt(reads, db)