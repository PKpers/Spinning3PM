#!/usr/bin/env python3


import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from fractions import Fraction
import sys

### Functions
def format_number(x):
    """Convert float to clean rational string"""
    frac = Fraction(float(x)).limit_denominator()
    if frac.denominator == 1:
        return str(frac.numerator)
    else:
        return f"{frac.numerator}/{frac.denominator}"
#
def normalize_region(v):
    """
    Normalize region vector v = (λ, x1, ..., xn)
    using Smirnov-Pak equivalence: v ~ v + (0,A,...,A).
    """
    lam = v[0]
    rest = np.array(v[1:], dtype=int)
    shift = rest.min()
    new_rest = rest - shift
    return tuple([lam] + new_rest.tolist())
#


### Main

# Get filename from command line argument
if len(sys.argv) != 2:
    print("Usage: python script.py <filename>")
    sys.exit(1)

filename = sys.argv[1]




# Read points from file
try:
    points = np.loadtxt(filename,dtype=int, delimiter=",")
    print(f"Successfully loaded {len(points)} points from {filename}...")
    print("constructing the convex hull...")
    hull = ConvexHull(points)
except FileNotFoundError:
    print(f"Error: File '{filename}' not found")
    sys.exit(1)
except Exception as e:
    print(f"Error reading file: {e}")
    sys.exit(1)

bottom_facets = []
inward_normals = []

print("determining the botom facets of the convex hull ... ")
for simplex, eq in zip(hull.simplices, hull.equations):
    a = eq[:-1]   # outward normal
    b = eq[-1]

    if a[0] < 0:  # outward λ-component negative => bottom facet
        bottom_facets.append(simplex)

        # inward normal
        n_in = -a
#
        # normalize so that λ=1
        n_in /= n_in[0]

        inward_normals.append(n_in)
print('done !')

#normalized = sorted(set([normalize_region(v) for v in inward_normals]))

# --- Plot hull with bottom facets ---
#plt.figure(figsize=(6,6))
#plt.plot(points[:,0], points[:,1], 'o')

#for simplex in hull.simplices:
#    plt.plot(points[simplex,0], points[simplex,1], 'k--', alpha=0.5)

#for simplex in bottom_facets:
#    plt.plot(points[simplex,0], points[simplex,1], 'r-', lw=3)

# show inward normals at facet midpoints
#for simplex, n in zip(bottom_facets, inward_normals):
#    midpoint = points[simplex].mean(axis=0)
#    plt.arrow(midpoint[0], midpoint[1], n[0]*0.1, n[1]*0.1,
#              head_width=0.02, color="blue")

#plt.xlabel("λ")
#plt.ylabel("x1")
#plt.title("Bottom facets with inward normals (λ=1)")
#plt.show()

# Process normals
normals_clean = []
for n in inward_normals:
    clean_tuple = tuple(format_number(component) for component in n)
    normals_clean.append(clean_tuple)

# Remove duplicates and sort
normals_clean = sorted(set(normals_clean))

print("found {} regions".format(len(normals_clean)))
for n in normals_clean:
    print(n)


'''
normals = [tuple(n) for n in inward_normals]

# deduplicate
normals = sorted(set(normals))

print("Inward normals (λ=1):")
#new = [tuple(row) for row in inward_normals]
#uniques = np.unique(new,axis=0)
for n in normals:
    print(n)
'''


