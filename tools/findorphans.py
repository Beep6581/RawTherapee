#!/usr/bin/env python3

import clang.cindex
import subprocess
import sys

index = clang.cindex.Index.create()
procevents = index.parse('rtengine/procevents.h',args=['-x', 'c++'])

if(1):
    for chld in procevents.cursor.get_children():
        if(chld.displayname == 'rtengine'):
            for c in chld.get_children():
                if(c.displayname == 'ProcEventCode'):
                    for pec in c.get_children():
                        #print(pec.kind, pec.displayname, pec.enum_value)
                        #print(pec.displayname, file=sys.stderr)
                        grp1 = subprocess.Popen(('grep', '-ro', '--exclude=procevents.h', '--exclude-dir=.git', pec.displayname), stdout=subprocess.PIPE)
                        wcr1 = subprocess.check_output(('wc', '-l'), stdin=grp1.stdout)
                        grp1.wait()
                        grp2 = subprocess.Popen(('grep', '-ro', '--exclude=procevents.h', '--exclude=refreshmap.cc', '--exclude-dir=.git', pec.displayname), stdout=subprocess.PIPE)
                        wcr2 = subprocess.check_output(('wc', '-l'), stdin=grp2.stdout)
                        grp2.wait()
                        print(pec.enum_value, pec.displayname,int(wcr1), int(wcr2))

with open('rtdata/languages/default', 'r') as deflang:
    for line in deflang:
        if(line[0] == '#'):
            continue
        if(line[0:2] == '//'):
            continue
        if(line[0:2] == '/*'):
            #our language files support comment blocks?????????????????????????????
            #or is this commented block bogus?
            continue
        if(line[0:2] == '*/'):
            continue
        else:
            stringid = line.split(';')[0]
            if(stringid.startswith('HISTORY_MSG')):
               continue
            #print(stringid, file=sys.stderr)
            grp1 = subprocess.Popen(('grep', '-ro', '--exclude-dir=languages', '--exclude-dir=.git', stringid), stdout=subprocess.PIPE)
            wcr1 = subprocess.check_output(('wc', '-l'), stdin=grp1.stdout)
            grp1.wait()
            print(stringid, int(wcr1))
