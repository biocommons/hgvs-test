# -*- coding: utf-8 -*-
from __future__ import absolute_import, division, print_function, unicode_literals

import recordtype
from hgvs.exceptions import HGVSError


class PosEdit(recordtype.recordtype('PosEdit', [('pos', None), ('edit', None), ('uncertain', False)])):
    """
    represents a **simple** variant, consisting of a single position and edit pair
    """

    def __str__(self):
        rv = str(self.edit) if self.pos is None else '{self.pos}{self.edit}'.format(self=self)
        if self.uncertain:
            if self.edit in ['0', '']:
                rv = rv + '?'
            else:
                rv = '(' + rv + ')'
        return rv

    def _set_uncertain(self):
        """sets the uncertain flag to True; used primarily by the HGVS grammar

        :returns: self
        """
        self.uncertain = True
        return self


class PosEdits:
    """
    represents a list of PosEdit for ComplexVariant
    """
    
    def __init__(self, posedits=[], uncertain=False):
        if type(posedits) is not list:
            raise HGVSError('posedits must be a list')
        self.posedits = posedits
        self.uncertain = uncertain
    
    def __getitem__(self, i):
        return self.posedits[i]

    def __setitem__(self, i, value):
        self.posedits[i] = value
    
    def __iter__(self):
        for posedit in self.posedits:
            yield posedit

    def __len__(self):
        return len(self.posedits)
    
    def __add__(self, set):
        return PosEdits(self.posedits + set.posedits, self.uncertain)
    
    @property
    def pos(self):
        return [posedit.pos for posedit in self.posedits]
    
    @property
    def edit(self):
        return [posedit.edit for posedit in self.posedits]
    
    @property
    def uncertain(self):
        return [posedit.uncertain for posedit in self.posedits]
    
    @property
    def phase_uncertain(self):
        return self.uncertain


## <LICENSE>
## Copyright 2014 HGVS Contributors (https://bitbucket.org/biocommons/hgvs)
## 
## Licensed under the Apache License, Version 2.0 (the "License");
## you may not use this file except in compliance with the License.
## You may obtain a copy of the License at
## 
##     http://www.apache.org/licenses/LICENSE-2.0
## 
## Unless required by applicable law or agreed to in writing, software
## distributed under the License is distributed on an "AS IS" BASIS,
## WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
## See the License for the specific language governing permissions and
## limitations under the License.
## </LICENSE>