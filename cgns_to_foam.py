#!/usr/bin/env python3
"""
cgns_to_foam.py  —  Convert WAND CGNS (unstructured HEXA_8, PointList BCs)
                    to OpenFOAM polyMesh.

Usage:
    python3 cgns_to_foam.py  <file.cgns>  [--scale 0.0254]  [--output constant/polyMesh]

Scale: WAND outputs in inches by default; 0.0254 converts to metres.
"""

import sys, os, argparse
import numpy as np
import h5py

# CGNS element type IDs
HEXA_8 = 17
QUAD_4 = 7


def dset(grp):
    """Read the ' data' dataset from a CGNS node group."""
    return grp[' data'][()]


def convert(cgns_path, output_dir, scale):
    print(f"Reading {cgns_path}")
    os.makedirs(output_dir, exist_ok=True)

    with h5py.File(cgns_path, 'r') as f:
        # Find CGNSBase (first Group in root that contains zone groups, not CGNSLibraryVersion)
        base = next(v for v in f.values()
                    if isinstance(v, h5py.Group) and v.name.split('/')[-1] != 'CGNSLibraryVersion'
                    and any(isinstance(c, h5py.Group) for c in v.values()))

        # Find the zone (first child Group that has GridCoordinates)
        zone = next(v for v in base.values()
                    if isinstance(v, h5py.Group) and 'GridCoordinates' in v)
        print(f"Zone: {zone.name.split('/')[-1]}")

        # ----------------------------------------------------------------
        # 1. Points
        # ----------------------------------------------------------------
        gc = zone['GridCoordinates']
        X = dset(gc['CoordinateX']).astype(np.float64) * scale
        Y = dset(gc['CoordinateY']).astype(np.float64) * scale
        Z = dset(gc['CoordinateZ']).astype(np.float64) * scale
        n_points = len(X)
        print(f"  {n_points} nodes")

        # ----------------------------------------------------------------
        # 2. Cells — read MixedElements, keep only HEXA_8
        # ----------------------------------------------------------------
        # ElementConnectivity stores: [elem_type, n0..n7, elem_type, n0..n7, ...]
        # ElementStartOffset[i] = start of element i in ElementConnectivity
        eg = zone['MixedElements']
        conn   = dset(eg['ElementConnectivity']).astype(np.int64)  # 1-indexed node IDs
        offset = dset(eg['ElementStartOffset']).astype(np.int64)

        n_cells = len(offset) - 1
        cells = np.zeros((n_cells, 8), dtype=np.int64)
        for ci in range(n_cells):
            s = offset[ci]
            etype = conn[s]
            if etype != HEXA_8:
                raise ValueError(f"Cell {ci} has element type {etype}, expected HEXA_8 ({HEXA_8})")
            # nodes are 1-indexed in CGNS → convert to 0-indexed
            cells[ci] = conn[s+1:s+9] - 1

        print(f"  {n_cells} HEXA_8 cells")

        # ----------------------------------------------------------------
        # 3. Boundary patches from ZoneBC (PointList = node indices)
        # ----------------------------------------------------------------
        # We build a node→{node set} membership map per patch.
        # A face belongs to a patch if ALL 4 of its nodes are in that patch's node set.
        patch_names = []
        patch_node_sets = []

        if 'ZoneBC' in zone:
            for pname, pgrp in zone['ZoneBC'].items():
                if not isinstance(pgrp, h5py.Group):
                    continue
                if 'PointList' not in pgrp:
                    print(f"  WARNING: patch {pname} has no PointList — skipping")
                    continue
                nodes_1idx = dset(pgrp['PointList']).flatten()
                node_set = set((nodes_1idx - 1).tolist())   # 0-indexed
                clean = sanitise_name(pname, patch_names)
                patch_names.append(clean)
                patch_node_sets.append(node_set)
                print(f"  Patch '{pname}' → '{clean}': {len(node_set)} nodes")

        # ----------------------------------------------------------------
        # 4. Build hex faces and classify internal / boundary
        # ----------------------------------------------------------------
        # OpenFOAM hex face ordering (6 faces per hex, each a quad):
        # Faces with outward normals for standard hex (nodes 0-7):
        #   0: bottom  (0,3,2,1)  -K
        #   1: top     (4,5,6,7)  +K
        #   2: front   (0,1,5,4)  -J
        #   3: back    (3,7,6,2)  +J
        #   4: left    (0,4,7,3)  -I
        #   5: right   (1,2,6,5)  +I
        HEX_FACES = [
            (0, 3, 2, 1),
            (4, 5, 6, 7),
            (0, 1, 5, 4),
            (3, 7, 6, 2),
            (0, 4, 7, 3),
            (1, 2, 6, 5),
        ]

        # face_map: canonical_key (frozenset) → dict
        #   first visit:  {'owner': cell, 'nodes': [4 ints], 'neighbour': -1}
        #   second visit: neighbour set
        face_map = {}

        print("Building face map...")
        for ci, hex_pts in enumerate(cells):
            for face_local in HEX_FACES:
                fn = [int(hex_pts[li]) for li in face_local]
                key = frozenset(fn)
                if key not in face_map:
                    face_map[key] = {'owner': ci, 'neighbour': -1, 'nodes': fn}
                else:
                    entry = face_map[key]
                    if entry['neighbour'] != -1:
                        print(f"WARNING: face {key} shared by >2 cells")
                        continue
                    nb = ci
                    if nb < entry['owner']:
                        # make smaller index the owner; flip normal
                        entry['neighbour'] = entry['owner']
                        entry['owner'] = nb
                        entry['nodes'] = entry['nodes'][::-1]
                    else:
                        entry['neighbour'] = nb

        print(f"  {len(face_map)} total unique faces")

        # Classify faces
        # A face is a boundary face if it has no neighbour (neighbour == -1).
        # Among boundary faces, assign to a patch if all 4 nodes ∈ patch node set.

        internal_faces  = []          # (owner, neighbour, nodes)
        patch_faces     = [[] for _ in patch_names]   # (owner, nodes) per patch
        unassigned      = []

        for key, entry in face_map.items():
            if entry['neighbour'] != -1:
                internal_faces.append((entry['owner'], entry['neighbour'], entry['nodes']))
            else:
                fn_set = key
                assigned = False
                for pi, nset in enumerate(patch_node_sets):
                    if fn_set.issubset(nset):
                        patch_faces[pi].append((entry['owner'], entry['nodes']))
                        assigned = True
                        break
                if not assigned:
                    unassigned.append(entry)

        internal_faces.sort(key=lambda x: (x[0], x[1]))

        print(f"  {len(internal_faces)} internal faces")
        for pi, pf in enumerate(patch_faces):
            print(f"  Patch '{patch_names[pi]}': {len(pf)} faces")
        if unassigned:
            print(f"  WARNING: {len(unassigned)} boundary faces not assigned to any patch")

        # ----------------------------------------------------------------
        # 5. Assemble face / owner / neighbour lists
        # ----------------------------------------------------------------
        all_faces = []
        all_owner = []
        all_neigh = []

        for owner, neighbour, nodes in internal_faces:
            all_faces.append(nodes)
            all_owner.append(owner)
            all_neigh.append(neighbour)

        n_internal = len(all_faces)
        patch_starts = []

        for pi, pf in enumerate(patch_faces):
            patch_starts.append(len(all_faces))
            for owner, nodes in pf:
                all_faces.append(nodes)
                all_owner.append(owner)

        # ----------------------------------------------------------------
        # 6. Write polyMesh
        # ----------------------------------------------------------------
        hdr = """\
/*--------------------------------*- C++ -*----------------------------------*\\
  =========                 |
  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\\\    /   O peration     | Website:  https://openfoam.org
    \\\\  /    A nd           | Version:  13
     \\\\/     M anipulation  |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    format      ascii;
    class       {cls};
    location    "constant/polyMesh";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""

        # points
        path = os.path.join(output_dir, 'points')
        print(f"Writing {path}")
        with open(path, 'w') as fh:
            fh.write(hdr.format(cls='vectorField', obj='points'))
            fh.write(f"{n_points}\n(\n")
            for x, y, z in zip(X, Y, Z):
                fh.write(f"({x:.10g} {y:.10g} {z:.10g})\n")
            fh.write(")\n")

        # faces
        path = os.path.join(output_dir, 'faces')
        print(f"Writing {path}")
        with open(path, 'w') as fh:
            fh.write(hdr.format(cls='faceList', obj='faces'))
            fh.write(f"{len(all_faces)}\n(\n")
            for nodes in all_faces:
                fh.write(f"4({nodes[0]} {nodes[1]} {nodes[2]} {nodes[3]})\n")
            fh.write(")\n")

        # owner
        path = os.path.join(output_dir, 'owner')
        print(f"Writing {path}")
        with open(path, 'w') as fh:
            fh.write(hdr.format(cls='labelList', obj='owner'))
            fh.write(f"{len(all_owner)}\n(\n")
            for o in all_owner:
                fh.write(f"{o}\n")
            fh.write(")\n")

        # neighbour (internal faces only)
        path = os.path.join(output_dir, 'neighbour')
        print(f"Writing {path}")
        with open(path, 'w') as fh:
            fh.write(hdr.format(cls='labelList', obj='neighbour'))
            fh.write(f"{n_internal}\n(\n")
            for n in all_neigh:
                fh.write(f"{n}\n")
            fh.write(")\n")

        # boundary
        path = os.path.join(output_dir, 'boundary')
        print(f"Writing {path}")
        with open(path, 'w') as fh:
            fh.write(hdr.format(cls='polyBoundaryMesh', obj='boundary'))
            n_patches = len(patch_names)
            fh.write(f"{n_patches}\n(\n")
            for pi, pname in enumerate(patch_names):
                n_faces = len(patch_faces[pi])
                start   = patch_starts[pi]
                ptype   = infer_type(pname)
                fh.write(f"    {pname}\n    {{\n")
                fh.write(f"        type            {ptype};\n")
                if ptype == 'cyclic':
                    nb_patch = cyclic_neighbour(pname, patch_names)
                    fh.write(f"        neighbourPatch  {nb_patch};\n")
                fh.write(f"        nFaces          {n_faces};\n")
                fh.write(f"        startFace       {start};\n")
                fh.write("    }\n")
            fh.write(")\n")

        print(f"\nDone → {output_dir}/")
        print(f"  {n_points} points | {n_cells} cells | "
              f"{n_internal} internal faces | "
              f"{len(all_faces)-n_internal} boundary faces")


def sanitise_name(name, existing):
    """Strip leading digits+underscore prefix; ensure unique valid OF word."""
    import re
    clean = re.sub(r'^[0-9]+_', '', name)
    if not clean or clean[0].isdigit():
        clean = 'patch_' + clean
    # deduplicate
    if clean in existing:
        i = 2
        while f"{clean}_{i}" in existing:
            i += 1
        clean = f"{clean}_{i}"
    return clean


def infer_type(name):
    n = name.upper()
    if 'PERIODIC' in n:
        return 'cyclic'
    if 'WALL' in n or 'SLIP' in n:
        return 'wall'
    return 'patch'


def cyclic_neighbour(name, all_names):
    n = name.upper()
    if 'SHADOW' in n:
        for other in all_names:
            if other != name and 'PERIODIC' in other.upper() and 'SHADOW' not in other.upper():
                return other
    else:
        for other in all_names:
            if other != name and 'PERIODIC' in other.upper() and 'SHADOW' in other.upper():
                return other
    return 'UNKNOWN'


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cgns', help='Input CGNS file')
    ap.add_argument('--scale', type=float, default=0.0254,
                    help='Coordinate scale factor (default 0.0254 = inches→metres)')
    ap.add_argument('--output', default='constant/polyMesh',
                    help='Output directory (default: constant/polyMesh)')
    args = ap.parse_args()
    convert(args.cgns, args.output, args.scale)


if __name__ == '__main__':
    main()
