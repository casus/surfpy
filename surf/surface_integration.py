"""
Created on Wed Dec  1 08:22:31 2023

Copyright 2023 by Gentian Zavalani.
"""
# Standard library imports:
import numpy as np
from scipy.linalg import norm
# Surfpy imports:
from transform import transform
from read_mesh_data import read_mesh_data
from SimpleImplicitSurfaceProjection import SimpleImplicitSurfaceProjection
from quad_weights import quadwts
from surfmesh import SurfceMesh


def integration(fun_handle, ls_function, grad_func, mesh, interp_deg):
    """
    Compute integration of a function over curved triangles.

    Args:
        fun_handle: callable, optional
            Function to be evaluated on each quadrature point.
        ls_function: callable
            Zero-levelset function.
        grad_func: callable
            Gradient of zero-levelset function.
        mesh: str
            The file path to the MAT file containing mesh data.
        interp_deg: int
            Interpolation degree.

    Returns:
        Integration values for each curved triangle.
    """
    # Read mesh data
    vertices, faces = read_mesh_data(mesh)  # Array of vertex coordinates and array of face connectivity.
    n_faces = faces.shape[0]  # Number of faces in the mesh
    n = interp_deg  # Interpolation degree
    pnts_quad = np.zeros((1, 3), dtype=np.float64)

    # Transform Chebyshev points from [-1,1]^2 to the reference simplex
    generating_points = transform(interp_deg, duffy_transform=True)
    quad_ps = np.array([[
        (1.0 - generating_points[row1, 0] - generating_points[row1, 1]),
        generating_points[row1, 0], generating_points[row1, 1]
    ] for row1 in range(generating_points.shape[0])])

    x = np.zeros((n_faces, interp_deg, interp_deg))
    y = np.zeros((n_faces, interp_deg, interp_deg))
    z = np.zeros((n_faces, interp_deg, interp_deg))

    # Loop over each face in the mesh
    for fun_id in range(n_faces):
        pnts_p = np.array([[0.0] * 3 for pid in range(generating_points.shape[0])])
        for q in range(quad_ps.shape[0]):
            pnts_quad = (
                quad_ps[q, 0] * vertices[faces[fun_id, 0]]
                + quad_ps[q, 1] * vertices[faces[fun_id, 1]]
                + quad_ps[q, 2] * vertices[faces[fun_id, 2]]
            )

            pnts_p[q] = SimpleImplicitSurfaceProjection(ls_function, grad_func, pnts_quad)

        # Reshape and store the results
        x[fun_id] = pnts_p[:, 0].reshape(interp_deg, interp_deg)
        y[fun_id] = pnts_p[:, 1].reshape(interp_deg, interp_deg)
        z[fun_id] = pnts_p[:, 2].reshape(interp_deg, interp_deg)

    # Create a surface mesh from the computed points
    dom = SurfceMesh(x, y, z)
    n = dom.x[0].shape[1] 

    # Flatten and concatenate the mesh data for function evaluation
    x_l = np.concatenate([subarray.flatten() for subarray in dom.x]).reshape(-1, 1)
    y_l = np.concatenate([subarray.flatten() for subarray in dom.y]).reshape(-1, 1)
    z_l = np.concatenate([subarray.flatten() for subarray in dom.z]).reshape(-1, 1)

    # Evaluate the function on the mesh
    f_values = fun_handle(x_l, y_l, z_l).reshape(len(dom.x), n, n)

    # Compute weights and Jacobian square root
    wu = quadwts(n).reshape(n, 1)
    wv = quadwts(n)
    J_sqrt = np.sqrt(np.abs(dom.J))

    # Calculate the final integration result
    Σ = np.sum(f_values * np.kron(wv, wu) * J_sqrt)

    return Σ
