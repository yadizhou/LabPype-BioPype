# -*- coding: utf-8 -*-

from .widget import *

ANCHORS = [
    (0, ANCHOR_SEQ_ANY, ANCHOR_SEQ_ANY),
    (0, ANCHOR_SEQ_DNA, ANCHOR_SEQ_DNA),
    (0, ANCHOR_SEQ_RNA, ANCHOR_SEQ_RNA),
    (0, ANCHOR_SEQ_PROTEIN, ANCHOR_SEQ_PROTEIN),

    (0, ANCHOR_SEQ_ANY, ANCHOR_SEQ_ANY_LIST),
    (0, ANCHOR_SEQ_DNA, ANCHOR_SEQ_DNA_LIST),
    (0, ANCHOR_SEQ_RNA, ANCHOR_SEQ_RNA_LIST),
    (0, ANCHOR_SEQ_PROTEIN, ANCHOR_SEQ_PROTEIN_LIST),

    (0, ANCHOR_REC_ANY, ANCHOR_REC_ANY),
    (0, ANCHOR_REC_DNA, ANCHOR_REC_DNA),
    (0, ANCHOR_REC_RNA, ANCHOR_REC_RNA),
    (0, ANCHOR_REC_PROTEIN, ANCHOR_REC_PROTEIN),

    (0, ANCHOR_REC_ANY, ANCHOR_REC_ANY_LIST),
    (0, ANCHOR_REC_DNA, ANCHOR_REC_DNA_LIST),
    (0, ANCHOR_REC_RNA, ANCHOR_REC_RNA_LIST),
    (0, ANCHOR_REC_PROTEIN, ANCHOR_REC_PROTEIN_LIST),

    (0, ANCHOR_REC_ANY, ANCHOR_SEQ_LIKE),
    (0, ANCHOR_SEQ_ANY, ANCHOR_SEQ_LIKE),
    (0, ANCHOR_REC_ANY_LIST, ANCHOR_SEQ_LIKE),
    (0, ANCHOR_SEQ_ANY_LIST, ANCHOR_SEQ_LIKE),
]

WIDGETS = [
    "File",
    ("#aaff55", LoadFile, "icon/Load.png"),
    ("#55ffaa", SaveFile, "icon/Save.png"),
    "Input",
    ("#ffbf00", Record, "icon/Sequence.png"),
    "Viewer",
    ("#ffc0ff", SequenceViewer, "icon/Viewer.png"),
]