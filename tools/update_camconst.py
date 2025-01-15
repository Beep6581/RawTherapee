#!/usr/bin/env python

from __future__ import print_function
import json
import argparse
import subprocess
import re
import os
import glob
import time


def getopts():
    p = argparse.ArgumentParser()
    p.add_argument('-c', '--cammatrices', required=True)
    p.add_argument('dir')
    p.add_argument('--dcamprof-path', default='.')
    p.add_argument('-v', '--version')
    return p.parse_args()


class Camconst(object):
    def __init__(self, comment_lines):
        self.comment = comment_lines
        self.data = {}
        self.aliases = {}

    def __getitem__(self, key):
        return self.data[key]

    def __contains__(self, key):
        return key in self.data

    def __setitem__(self, key, value):
        self.data[key] = value

# end of class Camconst


def load_camconst(fname):
    jsondata = []
    header = []
    with open(fname) as f:
        ok = True
        for line in f:
            l = line.strip()
            if l.startswith('/*'):
                ok = False
                header.append(l[2:].lstrip())
            elif l.endswith('*/'):
                header.append(l[:-2].rstrip())
                ok = True
            elif ok:
                jsondata.append(re.sub('//.*$', '', line))
            else:
                header.append(l)
    d = json.loads("\n".join(jsondata))
    res = Camconst(header)
    for entry in d['camera_constants']:
        try:
            camera = entry['make_model']
            if isinstance(camera, list):
                camera, rest = camera[0], camera
                res.aliases[camera] = rest
            res[camera] = entry['dcraw_matrix']
        except KeyError:
            pass
    return res


def dump_camconst(outname, camconst):
    with open(outname, 'w') as out:
        indent = ' ' * 4
        pr = out.write
        if camconst.comment:
            pr('/*\n')
            for line in camconst.comment:
                pr(line)
                pr('\n')
            pr('*/\n')
        pr('{"camera_constants": [\n')
        keysep = indent
        for key in sorted(camconst.data):
            pr(keysep)
            camera = key
            if key in camconst.aliases:
                camera = camconst.aliases[key]
            pr('{\n%s%s"make_model" : %s' % (indent, indent,
                                             json.dumps(camera)))
            sep = ',\n%s%s' % (indent, indent)
            pr('%s"dcraw_matrix" : %s' % (sep, json.dumps(camconst.data[key])))
            pr('\n%s}' % indent)
            keysep = ',\n%s' % indent
        pr('\n]}\n')


def extract_matrix(camconst, opts, filename):
    p = subprocess.Popen([os.path.join(opts.dcamprof_path, 'dcamprof'),
                          'dcp2json', filename], stdout=subprocess.PIPE)
    out, err = p.communicate()
    try:
        profile = json.loads(out)
        camera = profile['UniqueCameraModel'].upper()
        if camera in camconst:
            return None
        matrix = None
        for i in '1', '2':
            ill = 'CalibrationIlluminant' + i
            m = 'ColorMatrix' + i
            if ill in profile and m in profile and profile[ill] == 'D65':
                matrix = profile[m]
                break
        if matrix is not None and len(matrix) == 3 and len(matrix[0]) == 3:
            return camera, sum(([int(e * 10000) for e in row]
                                for row in matrix), [])
    except ValueError:
        return None


def main():
    opts = getopts()
    cammatrices = load_camconst(opts.cammatrices)
    updated = []
    for dcp in glob.glob(os.path.join(opts.dir, '*.dcp')):
        res = extract_matrix(cammatrices, opts, dcp)
        if res is not None:
            camera, matrix = res
            cammatrices[camera] = matrix
            print('Updated matrix for %s: %s' % (camera, matrix))
            updated.append(camera)
    if updated:
        info = ''
        if opts.version:
            info = ' (with Adobe DNG converter %s)' % opts.version
        cammatrices.comment.append('Updated on %s%s:' % (time.asctime(), info))
        for cam in sorted(updated):
            cammatrices.comment.append('  %s' % cam)
    dump_camconst(opts.cammatrices, cammatrices)


if __name__ == '__main__':
    main()
