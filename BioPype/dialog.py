# -*- coding: utf-8 -*-

from labpype.widget import Dialog

import wx

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO


class SequenceViewer(Dialog):
    SIZE = wx.Size(540, 400)

    def Initialize(self, Sizer):
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        self.List = self.AddListBox(sizer, width=120, onClick=self.OnSelect)
        self.Seq = self.AddTextCtrl(sizer, style=wx.TE_READONLY, font="FIXED")
        Sizer.Add(sizer, 1, wx.EXPAND)
        self.seqs = []

    def OnSelect(self):
        if self.List.HasSelection():
            s = self.seqs[self.List.GetSelection()][1]
            self.Seq.SetValue(str(s.seq) if isinstance(s, SeqRecord) else str(s))

    def GetData(self):
        self.seqs = []
        self.List.SetSelection(-1)
        self.List.ReDraw()
        self.Seq.SetValue("")
        if self.Widget["ANY"] is not None:
            for i in self.Widget["ANY"]:
                if i is not None:
                    if isinstance(i, (tuple, list)):
                        for j in i:
                            self.seqs.append((j.name, j))
                    else:
                        self.seqs.append((i.name, i))
        self.List.SetData([(i[0],) for i in self.seqs])
