#!/usr/bin/env python3
"""Freeze explicitly (--freeze) or rebuild a release from hash-checked local sources.
Maintainer dependencies: Python, numpy/scipy; R with jsonlite. Never downloads.
Run from the grip repository: python data-raw/build-suitesparse.py DATA_ROOT [--freeze].
"""
import argparse, gzip, hashlib, io, json, subprocess, tarfile, tempfile
from pathlib import Path
import numpy as np
from scipy.io import mmread
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

PACKAGE = Path(__file__).resolve().parents[1]
DEST = PACKAGE / 'inst/extdata/suitesparse'
SELECTION = DEST / 'selection-30_sept_2026.json'
SIX = ['HB/dwt_66', 'Newman/lesmis', 'HB/dwt_307', 'HB/494_bus', 'HB/dwt_1005', 'HB/1138_bus']
def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()
def checked(path, expected):
    if sha(path) != expected:
        raise ValueError(f'SHA-256 mismatch: {path}')
    return path

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('root', type=Path)
    parser.add_argument('--freeze', action='store_true', help='Explicitly replace release selection from current gallery')
    args = parser.parse_args(); root = args.root.resolve(); DEST.mkdir(exist_ok=True)
    if args.freeze:
        manifest = json.loads((root/'viewer_manifest.json').read_text())
        selection = dict(schema_version=1, release='30_sept_2026', gallery_manifest_sha256=sha(root/'viewer_manifest.json'), graphs=[])
        for e in sorted(manifest['graphs'], key=lambda e: e['id']):
            path = checked(root/e['file']['path'], e['file']['sha256'])
            g = json.loads(path.read_text())
            assert g['graph_sha256'] == e['graph_sha256']
            archive = Path(g['archive_file'])
            checked(archive, g['archive_sha256'])
            record = dict(id=e['id'], gallery=True, original_six=e['id'] in SIX,
                          graph_file=e['file'], graph_sha256=g['graph_sha256'],
                          archive=dict(path=str(archive.relative_to(root)), sha256=g['archive_sha256']))
            c = g.get('component_extraction')
            if c:
                mapping = checked(Path(c['vertex_mapping_file']), c['vertex_mapping_sha256'])
                original = Path(c['original_graph_file'])
                record['mapping'] = dict(path=str(mapping.relative_to(root)), sha256=sha(mapping))
                record['original_graph'] = dict(path=str(original.relative_to(root)), sha256=sha(original))
            selection['graphs'].append(record)
        for id in SIX:
            p = PACKAGE/'inst/extdata/zheng-graphs'/f'{id.split("/")[1]}.mtx.gz'
            entry = next((e for e in selection['graphs'] if e['id']==id), None)
            if entry is None:
                entry = dict(id=id, gallery=False, original_six=True)
                selection['graphs'].append(entry)
            entry['original_matrix'] = dict(path=str(p.relative_to(PACKAGE)), sha256=sha(p))
        labels = PACKAGE/'inst/extdata/zheng-graphs/lesmis_nodename.txt'
        selection['labels'] = dict(path=str(labels.relative_to(PACKAGE)), sha256=sha(labels))
        selection['graphs'].sort(key=lambda e: e['id'])
    else:
        selection = json.loads(SELECTION.read_text())
    checked(PACKAGE/selection['labels']['path'], selection['labels']['sha256'])
    output = []
    for e in selection['graphs']:
        id=e['id']; g=None; component=None
        if 'original_matrix' in e:
            checked(PACKAGE/e['original_matrix']['path'],e['original_matrix']['sha256'])
        if e['gallery']:
            g=json.loads(checked(root/e['graph_file']['path'], e['graph_file']['sha256']).read_text())
            assert g['graph_id']==id and g['graph_sha256']==e['graph_sha256']
            identity={k:g[k] for k in ('recipe','vertex_ids','edges','edge_length')}
            assert hashlib.sha256(json.dumps(identity,sort_keys=True,allow_nan=False,separators=(',', ':')).encode()).hexdigest()==e['graph_sha256']
            component=g.get('component_extraction')
            source_id=g['source']['graph_id']
            with tarfile.open(checked(root/e['archive']['path'],e['archive']['sha256'])) as t:
                members=[m for m in t.getmembers() if m.isfile() and Path(m.name).name==source_id.split('/')[1]+'.mtx']
                assert len(members)==1
                raw=t.extractfile(members[0]).read()
        else:
            source_id=id
            raw=gzip.decompress(checked(PACKAGE/e['original_matrix']['path'],e['original_matrix']['sha256']).read_bytes())
        matrix=mmread(io.BytesIO(raw), spmatrix=True).tocoo(); matrix.sum_duplicates(); matrix.eliminate_zeros()
        assert np.isfinite(matrix.data).all()
        rows,cols=matrix.shape; bipartite=rows!=cols; original_n=rows+cols if bipartite else rows
        pairs=sorted({(min(int(a),int(b)),max(int(a),int(b))) for a,b in zip(matrix.row,matrix.col+rows if bipartite else matrix.col) if bipartite or a!=b})
        original_indices=list(range(1,original_n+1))
        if component:
            mapping=json.loads(checked(root/e['mapping']['path'],e['mapping']['sha256']).read_text())
            original=json.loads(checked(root/e['original_graph']['path'],e['original_graph']['sha256']).read_text())
            assert original['graph_sha256']==component['original_graph_sha256']
            assert original['edges']==[list(p) for p in pairs]
            original_indices=component['original_vertex_indices_one_based']
            assert original_indices==sorted(set(original_indices))
            assert [v['original_vertex'] for v in mapping]==original_indices
            assert [v['vertex'] for v in mapping]==list(range(1,len(mapping)+1))
            assert [v['original_vertex_id'] for v in mapping]==[original['vertex_ids'][i-1] for i in original_indices]
            remap={old-1:new for new,old in enumerate(original_indices)}
            pairs=[(remap[a],remap[b]) for a,b in pairs if a in remap and b in remap]
            assert component['original_vertices']==original_n
        n=len(original_indices)
        edges=np.array(pairs,dtype=int).reshape(-1,2)
        adjacency=csr_matrix((np.ones(2*len(edges)),(np.r_[edges[:,0],edges[:,1]],np.r_[edges[:,1],edges[:,0]])),shape=(n,n))
        count,labels=connected_components(adjacency,directed=False)
        if g:
            assert g['edges']==edges.tolist() and g['n_vertices']==n and g['n_edges']==len(edges) and g['n_components']==count, id
            if component:
                assert count==1 and n/original_n>=.95
        header=[]
        for line in raw.decode().splitlines():
            if not line.startswith('%'): break
            header.append(line)
        metadata=dict(name=id.split('/')[1]+'_unweighted_graph', graph_id=id, original_graph_id=source_id, source_url='https://sparse.tamu.edu/'+source_id,
                      download_url='https://sparse.tamu.edu/MM/'+source_id+'.tar.gz',
                      license='CC BY 4.0', license_url='https://creativecommons.org/licenses/by/4.0/',
                      source_header=header, matrix_sha256=hashlib.sha256(raw).hexdigest(),
                      conversion_type='bipartite' if bipartite else ('largest_component' if component else 'square'),
                      conversion_recipe='Sum duplicate entries; discard numerical zeros. '+('Rows and columns form separate vertex sets.' if bipartite else 'Remove diagonal; take union of off-diagonal support with its transpose.')+(' Induce largest component, retaining source vertex order.' if component else ' Retain all source vertices.')+' Store each undirected edge once; use unit traversal lengths.',
                      original_vertices=original_n, retained_vertices=n, retained_fraction=n/original_n,
                      source_rows=rows, source_columns=cols, selection_release=selection['release'])
        if g:
            metadata.update(archive_sha256=e['archive']['sha256'],graph_json_sha256=e['graph_file']['sha256'],graph_sha256=e['graph_sha256'],source_recipe=g['recipe'])
        if component:
            metadata.update(original_graph_sha256=component['original_graph_sha256'],vertex_mapping_sha256=e['mapping']['sha256'])
        vertex_ids=g['vertex_ids'] if g else [f'v:{i+1}' for i in range(n)]
        assert len(vertex_ids)==n and len(set(vertex_ids))==n
        source_vertex_ids=[v['original_vertex_id'] for v in mapping] if component else vertex_ids
        vertex_labels=[str(i) for i in original_indices]
        partitions=['row']*n; source_indices=original_indices
        if bipartite:
            partitions=['row']*rows+['column']*cols
            source_indices=list(range(1,rows+1))+list(range(1,cols+1))
            vertex_labels=[f'{p}:{i}' for p,i in zip(partitions,source_indices)]
        if id=='Newman/lesmis': vertex_labels=(PACKAGE/selection['labels']['path']).read_text().splitlines()
        output.append(dict(id=id,n=n,edges=(edges+1).tolist(),vertex_labels=vertex_labels,
                           source_indices=source_indices,original_indices=original_indices,partitions=partitions,
                           vertex_ids=vertex_ids,source_vertex_ids=source_vertex_ids,
                           components=(labels+1).tolist(),n_components=int(count),provenance=metadata,
                           gallery=e['gallery'],original_six=e['original_six']))
    with tempfile.NamedTemporaryFile(mode='w',suffix='.json') as tmp:
        json.dump(output,tmp); tmp.flush()
        subprocess.run(['Rscript',str(PACKAGE/'data-raw/suitesparse_graphs.R'),tmp.name],check=True,cwd=PACKAGE)
    if args.freeze:
        SELECTION.write_text(json.dumps(selection, indent=2)+'\n')
    print(f'Validated {len(output)} graphs against source matrices and frozen checksums.')
if __name__=='__main__': main()
