# -*- coding: utf-8 -*-

from labpype.widget import Widget, ANCHOR_REGULAR
from labpype.widget.field import *
from . import dialog as Dlg

from collections import OrderedDict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO


# =================================================== Anchor Types =================================================== #
# @formatter:off

class ANCHOR_SEQ_LIKE(ANCHOR_REGULAR): pass

class ANCHOR_SEQ_ANY(ANCHOR_REGULAR): pass
class ANCHOR_SEQ_DNA(ANCHOR_REGULAR): pass
class ANCHOR_SEQ_RNA(ANCHOR_REGULAR): pass
class ANCHOR_SEQ_PROTEIN(ANCHOR_REGULAR): pass

class ANCHOR_SEQ_ANY_LIST(ANCHOR_REGULAR): pass
class ANCHOR_SEQ_DNA_LIST(ANCHOR_REGULAR): pass
class ANCHOR_SEQ_RNA_LIST(ANCHOR_REGULAR): pass
class ANCHOR_SEQ_PROTEIN_LIST(ANCHOR_REGULAR): pass

class ANCHOR_REC_ANY(ANCHOR_REGULAR): pass
class ANCHOR_REC_DNA(ANCHOR_REGULAR): pass
class ANCHOR_REC_RNA(ANCHOR_REGULAR): pass
class ANCHOR_REC_PROTEIN(ANCHOR_REGULAR): pass

class ANCHOR_REC_ANY_LIST(ANCHOR_REGULAR): pass
class ANCHOR_REC_DNA_LIST(ANCHOR_REGULAR): pass
class ANCHOR_REC_RNA_LIST(ANCHOR_REGULAR): pass
class ANCHOR_REC_PROTEIN_LIST(ANCHOR_REGULAR): pass

# class ANCHOR_DNA_LIST(ANCHOR_REGULAR): pass
# class ANCHOR_RNA_LIST(ANCHOR_REGULAR): pass
# class ANCHOR_PROTEIN_LIST(ANCHOR_REGULAR): pass
# @formatter:on

# ====================================================== Widget ====================================================== #
MAP_ALPHABET = OrderedDict((
    ("DNA", IUPAC.unambiguous_dna),
    ("RNA", IUPAC.unambiguous_rna),
    ("Protein", IUPAC.protein),
    ("Ambiguous DNA", IUPAC.ambiguous_dna),
    ("Ambiguous RNA", IUPAC.ambiguous_rna),
    ("Extended Protein", IUPAC.extended_protein),
))


# ======================================================== #
class Sequence(Widget):
    NAME = "New Sequence"
    DIALOG = {"ORIENTATION": "V", "SIZE": (540, 400)}
    INTERNAL = LineField(key="NAME", label="", hint="Input Sequence Name Here"), \
               ChoiceField(key="TYPE", label="", choices=list(MAP_ALPHABET.keys()), width=-1), \
               TextField(key="SEQ", label="", font="FIXED")
    OUTGOING = ANCHOR_SEQ_ANY

    def Name(self):
        return self["NAME"]

    def Task(self):
        seq = Seq(self["SEQ"], alphabet=MAP_ALPHABET[self["TYPE"]])
        seq.name = self["NAME"]
        return seq


# ======================================================== #
class Record(Widget):
    NAME = "New Record"
    DIALOG = {"ORIENTATION": "V", "SIZE": (540, 400)}
    INTERNAL = LineField(key="NAME", label="", hint="Input Record Name Here"), \
               ChoiceField(key="TYPE", label="", choices=list(MAP_ALPHABET.keys()), width=-1), \
               TextField(key="SEQ", label="", font="FIXED")
    OUTGOING = ANCHOR_REC_ANY

    def Name(self):
        return self["NAME"]

    def Task(self):
        return SeqRecord(Seq(self["SEQ"], alphabet=MAP_ALPHABET[self["TYPE"]]), name=self["NAME"])


# ======================================================== #
class LoadFile(Widget):
    NAME = "Load File"
    THREAD = True
    DIALOG = {"ORIENTATION": "V", "SIZE": (320, -1)}
    INTERNAL = ChoiceField(key="TYPE", label="", choices=list(MAP_ALPHABET.keys()), width=-1), \
               FileField(key="FILE", label="", mode="L", wildcard="FASTA files (*.fasta)|*.fasta|GenBank files (*.gbk)|*.gbk")
    OUTGOING = ANCHOR_REC_ANY_LIST

    def Name(self):
        if self["OUT"] is not None:
            return "%s %s Record%s" % (len(self["OUT"]), self["TYPE"], "s" if len(self["OUT"]) > 1 else "")

    def Task(self):
        o = []
        ft = None
        if self["FILE"].endswith(".fasta"):
            ft = "fasta"
        elif self["FILE"].endswith(".gbk"):
            ft = "genbank"
        for index,i in enumerate(SeqIO.parse(self["FILE"], ft, alphabet=MAP_ALPHABET[self["TYPE"]])):
            o.append(i)
            self.Checkpoint("%s" % (index + 1))
        return o


class SaveFile(Widget):
    NAME = "Save File"
    THREAD = True
    DIALOG = {"ORIENTATION": "V", "SIZE": (320, -1)}
    INTERNAL = ChoiceField(key="TYPE", label="Save as ...", choices=("FASTA", "GenBank"), width=-1), \
               FileField(key="FILE", label="Path ...", mode="S", wildcard="FASTA files (*.fasta)|*.fasta|GenBank files (*.gbk)|*.gbk")
    INCOMING = ANCHOR_SEQ_LIKE, "ANY", True, "LTB"

    def Task(self):
        ft = None
        if self["TYPE"] == "FASTA":
            ft = "fasta"
        elif self["TYPE"] == "GenBank":
            ft = "genbank"
        with open(self["FILE"], "w") as fo:
            for s in self["ANY"]:
                self.Checkpoint()
                if isinstance(s, Seq):
                    SeqIO.write(SeqRecord(s, name=s.name), fo, ft)
                else:
                    SeqIO.write(s, fo, ft)
        return True


# ======================================================== #
class SequenceViewer(Widget):
    NAME = "View Sequence"
    DIALOG = Dlg.SequenceViewer
    INCOMING = ANCHOR_SEQ_LIKE, "ANY", True, "LTB"

    def Task(self):
        return True

    def IsIncomingAvailable(self):
        for a in self.Incoming:
            if a.connected:
                return True
        return False
