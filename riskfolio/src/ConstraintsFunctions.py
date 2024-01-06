""""""  #
"""
Copyright (c) 2020-2024, Dany Cajas
All rights reserved.
This work is licensed under BSD 3-Clause "New" or "Revised" License.
License available at https://github.com/dcajasn/Riskfolio-Lib/blob/master/LICENSE.txt
"""

import numpy as np
import pandas as pd
import networkx as nx
import scipy.cluster.hierarchy as hr
from scipy.spatial.distance import squareform
import riskfolio.src.AuxFunctions as af
import riskfolio.src.DBHT as db


__all__ = [
    "assets_constraints",
    "factors_constraints",
    "assets_views",
    "factors_views",
    "assets_clusters",
    "hrp_constraints",
    "risk_constraint",
    "connection_matrix",
    "centrality_vector",
    "clusters_matrix",
    "connected_assets",
    "related_assets",
]


def assets_constraints(constraints, asset_classes):
    r"""
    Create the linear constraints matrixes A and B of the constraint
    :math:`Aw \geq B`.

    Parameters
    ----------
    constraints : DataFrame of shape (n_constraints, n_fields)
        Constraints matrix, where n_constraints is the number of constraints
        and n_fields is the number of fields of constraints matrix, the fields
        are:

        - Disabled: (bool) indicates if the constraint is enable.
        - Type: (str) can be: 'Assets', 'Classes', 'All Assets', 'Each asset in a class' and 'All Classes'.
        - Set: (str) if Type is 'Classes', 'Each asset in a class' or 'All Classes'specified the name of the asset's classes set.
        - Position: (str) the name of the asset or asset class of the constraint.
        - Sign: (str) can be '>=' or '<='.
        - Weight: (scalar) is the maximum or minimum weight of the absolute constraint.
        - Type Relative: (str) can be: 'Assets' or 'Classes'.
        - Relative Set: (str) if Type Relative is 'Classes' specified the name of the set of asset classes.
        - Relative: (str) the name of the asset or asset class of the relative constraint.
        - Factor: (scalar) is the factor of the relative constraint.

    asset_classes : DataFrame of shape (n_assets, n_cols)
        Asset's classes matrix, where n_assets is the number of assets and
        n_cols is the number of columns of the matrix where the first column
        is the asset list and the next columns are the different asset's
        classes sets.

    Returns
    -------
    A : nd-array
        The matrix A of :math:`Aw \geq B`.

    B : nd-array
        The matrix B of :math:`Aw \geq B`.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------
    ::

        import riskfolio as rp

        asset_classes = {'Assets': ['FB', 'GOOGL', 'NTFX', 'BAC', 'WFC', 'TLT', 'SHV'],
                         'Class 1': ['Equity', 'Equity', 'Equity', 'Equity', 'Equity',
                                     'Fixed Income', 'Fixed Income'],
                         'Class 2': ['Technology', 'Technology', 'Technology',
                                     'Financial', 'Financial', 'Treasury', 'Treasury'],}

        asset_classes = pd.DataFrame(asset_classes)
        asset_classes = asset_classes.sort_values(by=['Assets'])

        constraints = {'Disabled': [False, False, False, False, False, False, False],
                       'Type': ['Classes', 'Classes', 'Assets', 'Assets', 'Classes',
                                'All Assets', 'Each asset in a class'],
                       'Set': ['Class 1', 'Class 1', '', '', 'Class 2', '', 'Class 1'],
                       'Position': ['Equity', 'Fixed Income', 'BAC', 'WFC', 'Financial',
                                    '', 'Equity'],
                       'Sign': ['<=', '<=', '<=', '<=', '>=', '>=', '>='],
                       'Weight': [0.6, 0.5, 0.1, '', '', 0.02, ''],
                       'Type Relative': ['', '', '', 'Assets', 'Classes', '', 'Assets'],
                       'Relative Set': ['', '', '', '', 'Class 1', '', ''],
                       'Relative': ['', '', '', 'FB', 'Fixed Income', '', 'TLT'],
                       'Factor': ['', '', '', 1.2, 0.5, '', 0.4]}

        constraints = pd.DataFrame(constraints)


    The constraint looks like this:

    .. image:: images/Constraints.png

    It is easier to construct the constraints in excel and then upload to a
    dataframe.

    To create the matrixes A and B we use the following command:

    ::

        A, B = rp.assets_constraints(constraints, asset_classes)


    The matrixes A and B looks like this (all constraints were converted to a linear
    constraint):

    .. image:: images/AxB.png

    """

    if not isinstance(constraints, pd.DataFrame) and not isinstance(
        asset_classes, pd.DataFrame
    ):
        raise ValueError("constraints and asset_classes must be DataFrames")

    if constraints.shape[1] != 10:
        raise ValueError("constraints must have ten columns")

    n = len(constraints)
    m = len(asset_classes)
    data = constraints.fillna("")
    data = data.values.tolist()
    assetslist = asset_classes.iloc[:, 0].values.tolist()

    A = []
    B = []
    for i in range(0, n):
        if data[i][0] == False:
            if data[i][1] == "Assets":
                item = assetslist.index(data[i][3])
                if data[i][4] == ">=":
                    d = 1
                elif data[i][4] == "<=":
                    d = -1
                if data[i][5] != "":
                    A1 = [0] * m
                    A1[item] = d
                    A.append(A1)
                    B.append([data[i][5] * d])
                else:
                    A1 = [0] * m
                    A1[item] = 1
                    if data[i][6] == "Assets":
                        item2 = assetslist.index(data[i][8])
                        A2 = [0] * m
                        A2[item2] = 1
                    elif data[i][6] == "Classes":
                        A2 = np.where(
                            asset_classes[data[i][7]].values == data[i][8], 1, 0
                        )
                    A1 = ((np.array(A1) + np.array(A2) * data[i][9] * -1) * d).tolist()
                    A.append(A1)
                    B.append([0])
            elif data[i][1] == "All Assets":
                item = len(assetslist)
                if data[i][4] == ">=":
                    d = 1
                elif data[i][4] == "<=":
                    d = -1
                if data[i][5] != "":
                    A1 = np.identity(item) * d
                    A1 = A1.tolist()
                    B1 = np.ones((item, 1)) * d * data[i][5]
                    for i in range(0, item):
                        A.append(A1[i])
                        B.append(B1.tolist()[0])
                else:
                    A1 = np.identity(item)
                    if data[i][6] == "Assets":
                        item2 = assetslist.index(data[i][8])
                        A2 = np.zeros((item, item - 1))
                        A2 = np.insert(A2, item2 - 1, 1, axis=1)
                    elif data[i][6] == "Classes":
                        A1 = np.identity(item)
                        A2 = np.where(
                            asset_classes[data[i][7]].values == data[i][8], 1, 0
                        )
                        A2 = np.ones((item, item)) * np.array(A2)
                    A1 = ((np.array(A1) + np.array(A2) * data[i][9] * -1) * d).tolist()
                    for i in range(0, item):
                        A.append(A1[i])
                        B.append([0])
            elif data[i][1] == "Classes":
                if data[i][4] == ">=":
                    d = 1
                elif data[i][4] == "<=":
                    d = -1
                if data[i][5] != "":
                    A1 = np.where(asset_classes[data[i][2]].values == data[i][3], 1, 0)
                    A1 = np.array(A1) * d
                    A1 = A1.tolist()
                    A.append(A1)
                    B.append([data[i][5] * d])
                else:
                    A1 = np.where(asset_classes[data[i][2]].values == data[i][3], 1, 0)
                    if data[i][6] == "Assets":
                        item2 = assetslist.index(data[i][8])
                        A2 = [0] * m
                        A2[item2] = 1
                    elif data[i][6] == "Classes":
                        A2 = np.where(
                            asset_classes[data[i][7]].values == data[i][8], 1, 0
                        )
                    A1 = ((np.array(A1) + np.array(A2) * data[i][9] * -1) * d).tolist()
                    A.append(A1)
                    B.append([0])
            elif data[i][1] == "Each asset in a class":
                if data[i][4] == ">=":
                    d = 1
                elif data[i][4] == "<=":
                    d = -1
                if data[i][5] != "":
                    A1 = np.where(asset_classes[data[i][2]].values == data[i][3], 1, 0)
                    l = 0
                    for k in A1:
                        if k == 1:
                            A3 = [0] * m
                            A3[l] = 1 * d
                            A.append(A3)
                            B.append([data[i][5] * d])
                        l = l + 1
                else:
                    A1 = np.where(asset_classes[data[i][2]].values == data[i][3], 1, 0)
                    l = 0
                    for k in A1:
                        if k == 1:
                            A3 = [0] * m
                            A3[l] = 1
                            if data[i][6] == "Assets":
                                item2 = assetslist.index(data[i][8])
                                A2 = [0] * m
                                A2[item2] = 1
                            elif data[i][6] == "Classes":
                                A2 = np.where(
                                    asset_classes[data[i][7]].values == data[i][8], 1, 0
                                )
                            A3 = (
                                (np.array(A3) + np.array(A2) * data[i][9] * -1) * d
                            ).tolist()
                            A.append(A3)
                            B.append([0])
                        l = l + 1
            elif data[i][1] == "All Classes":
                if data[i][4] == ">=":
                    d = 1
                elif data[i][4] == "<=":
                    d = -1
                if data[i][5] != "":
                    for k in np.unique(asset_classes[data[i][2]].values):
                        A1 = np.where(asset_classes[data[i][2]].values == k, 1, 0) * d
                        A1 = A1.tolist()
                        A.append(A1)
                        B.append([data[i][5] * d])
                else:
                    for k in np.unique(asset_classes[data[i][2]].values):
                        A1 = np.where(asset_classes[data[i][2]].values == k, 1, 0)
                        if data[i][6] == "Assets":
                            item2 = assetslist.index(data[i][8])
                            A2 = [0] * m
                            A2[item2] = 1
                        elif data[i][6] == "Classes":
                            A2 = np.where(
                                asset_classes[data[i][7]].values == data[i][8], 1, 0
                            )
                        A3 = (
                            (np.array(A1) + np.array(A2) * data[i][9] * -1) * d
                        ).tolist()
                        A.append(A3)
                        B.append([0])

    A = np.array(A, ndmin=2, dtype=float)
    B = np.array(B, ndmin=2, dtype=float)

    return A, B


def factors_constraints(constraints, loadings):
    r"""
    Create the factors constraints matrixes C and D of the constraint
    :math:`Cw \geq D`.

    Parameters
    ----------
    constraints : DataFrame of shape (n_constraints, n_fields)
        Constraints matrix, where n_constraints is the number of constraints
        and n_fields is the number of fields of constraints matrix, the fields
        are:

        - Disabled: (bool) indicates if the constraint is enable.
        - Factor: (str) the name of the factor of the constraint.
        - Sign: (str) can be '>=' or '<='.
        - Value: (scalar) is the maximum or minimum value of the factor.

    loadings : DataFrame of shape (n_assets, n_features)
        The loadings matrix.

    Returns
    -------
    C : nd-array
        The matrix C of :math:`Cw \geq D`.

    D : nd-array
        The matrix D of :math:`Cw \geq D`.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------
    ::

        loadings = {'const': [0.0004, 0.0002, 0.0000, 0.0006, 0.0001, 0.0003, -0.0003],
                    'MTUM': [0.1916, 1.0061, 0.8695, 1.9996, 0.0000, 0.0000, 0.0000],
                    'QUAL': [0.0000, 2.0129, 1.4301, 0.0000, 0.0000, 0.0000, 0.0000],
                    'SIZE': [0.0000, 0.0000, 0.0000, 0.4717, 0.0000, -0.1857, 0.0000],
                    'USMV': [-0.7838, -1.6439, -1.0176, -1.4407, 0.0055, 0.5781, 0.0000],
                    'VLUE': [1.4772, -0.7590, -0.4090, 0.0000, -0.0054, -0.4844, 0.9435]}

        loadings = pd.DataFrame(loadings)

        constraints = {'Disabled': [False, False, False],
                       'Factor': ['MTUM', 'USMV', 'VLUE'],
                       'Sign': ['<=', '<=', '>='],
                       'Value': [0.9, -1.2, 0.3],
                       'Relative Factor': ['USMV', '', '']}

        constraints = pd.DataFrame(constraints)


    The constraint looks like this:

    .. image:: images/Constraints2.png

    It is easier to construct the constraints in excel and then upload to a
    dataframe.

    To create the matrixes C and D we use the following command:

    ::

        C, D = rp.factors_constraints(constraints, loadings)


    The matrixes C and D looks like this (all constraints were converted to a linear
    constraint):

    .. image:: images/CxD.png

    """

    if not isinstance(constraints, pd.DataFrame) and not isinstance(
        loadings, pd.DataFrame
    ):
        raise ValueError("constraints and loadings must be DataFrames")

    if constraints.shape[1] != 5:
        raise ValueError("constraints must have five columns")

    n = len(constraints)
    data = constraints.fillna("")
    data = data.values.tolist()

    C = []
    D = []
    for i in range(0, n):
        if data[i][0] == False:
            if data[i][2] == ">=":
                d = 1
            elif data[i][2] == "<=":
                d = -1
            C1 = loadings[data[i][1]].values
            if data[i][4] != "":
                C2 = loadings[data[i][4]].values
                C1 = C1 - C2
            C.append(C1 * d)
            D.append([data[i][3] * d])

    C = np.array(C, ndmin=2, dtype=float)
    D = np.array(D, ndmin=2, dtype=float)

    return C, D


def assets_views(views, asset_classes):
    r"""
    Create the assets views matrixes P and Q of the views :math:`Pw = Q`.

    Parameters
    ----------
    views : DataFrame of shape (n_views, n_fields)
        Constraints matrix, where n_views is the number of views
        and n_fields is the number of fields of views matrix, the fields
        are:

        - Disabled: (bool) indicates if the constraint is enable.
        - Type: (str) can be: 'Assets' or 'Classes'.
        - Set: (str) if Type is 'Classes' specified the name of the set of asset classes.
        - Position: (str) the name of the asset or asset class of the view.
        - Sign: (str) can be '>=' or '<='.
        - Return: (scalar) is the return of the view.
        - Type Relative: (str) can be: 'Assets' or 'Classes'.
        - Relative Set: (str) if Type Relative is 'Classes' specified the name of the set of asset classes.
        - Relative: (str) the name of the asset or asset class of the relative view.

    asset_classes : DataFrame of shape (n_assets, n_cols)
        Asset's classes matrix, where n_assets is the number of assets and
        n_cols is the number of columns of the matrix where the first column
        is the asset list and the next columns are the different asset's
        classes sets.

    Returns
    -------
    P : nd-array
        The matrix P that shows the relation among assets in each view.

    Q : nd-array
        The matrix Q that shows the expected return of each view.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------
    ::

        asset_classes = {'Assets': ['FB', 'GOOGL', 'NTFX', 'BAC', 'WFC', 'TLT', 'SHV'],
                         'Class 1': ['Equity', 'Equity', 'Equity', 'Equity', 'Equity',
                                      'Fixed Income', 'Fixed Income'],
                         'Class 2': ['Technology', 'Technology', 'Technology',
                                      'Financial', 'Financial', 'Treasury', 'Treasury'],}

        asset_classes = pd.DataFrame(asset_classes)
        asset_classes = asset_classes.sort_values(by=['Assets'])

        views = {'Disabled': [False, False, False, False],
                 'Type': ['Assets', 'Classes', 'Classes', 'Assets'],
                 'Set': ['', 'Class 2','Class 1', ''],
                 'Position': ['WFC', 'Financial', 'Equity', 'FB'],
                 'Sign': ['<=', '>=', '>=', '>='],
                 'Return': [ 0.3, 0.1, 0.05, 0.03 ],
                 'Type Relative': [ 'Assets', 'Classes', 'Assets', ''],
                 'Relative Set': [ '', 'Class 1', '', ''],
                 'Relative': ['FB', 'Fixed Income', 'TLT', '']}

        views = pd.DataFrame(views)


    The constraint looks like this:

    .. image:: images/Views.png

    It is easier to construct the constraints in excel and then upload to a
    dataframe.

    To create the matrixes P and Q we use the following command:

    ::

        P, Q = rp.assets_views(views, asset_classes)


    The matrixes P and Q looks like this:

    .. image:: images/PxQ.png

    """

    if not isinstance(views, pd.DataFrame) and not isinstance(
        asset_classes, pd.DataFrame
    ):
        raise ValueError("constraints and asset_classes must be DataFrames")

    if views.shape[1] != 9:
        raise ValueError("constraints must have nine columns")

    n = len(views)
    m = len(asset_classes)
    data = views.fillna("")
    data = data.values.tolist()
    assetslist = asset_classes.iloc[:, 0].values.tolist()

    P = []
    Q = []
    for i in range(0, n):
        valid = False
        if data[i][0] == False:
            if data[i][1] == "Assets":
                item = assetslist.index(data[i][3])
                if data[i][4] == ">=":
                    d = 1
                elif data[i][4] == "<=":
                    d = -1
                if data[i][5] != "":
                    P1 = [0] * m
                    P1[item] = 1
                    if data[i][6] == "Assets" and data[i][8] != "":
                        item2 = assetslist.index(data[i][8])
                        P2 = [0] * m
                        P2[item2] = 1
                        valid = True
                    elif (
                        data[i][6] == "Classes"
                        and data[i][7] != ""
                        and data[i][8] != ""
                    ):
                        P2 = np.where(
                            asset_classes[data[i][7]].values == data[i][8], 1, 0
                        )
                        P2 = P2 / np.sum(P2)
                        valid = True
                    elif data[i][6] == "" and data[i][7] == "" and data[i][8] == "":
                        P2 = [0] * m
                        valid = True
                    if valid == True:
                        P1 = ((np.array(P1) - np.array(P2)) * d).tolist()
                        P.append(P1)
                        Q.append([data[i][5] * d])
            elif data[i][1] == "Classes":
                if data[i][4] == ">=":
                    d = 1
                else:
                    d = -1
                if data[i][5] != "":
                    P1 = np.where(asset_classes[data[i][2]].values == data[i][3], 1, 0)
                    P1 = P1 / np.sum(P1)
                    if data[i][6] == "Assets" and data[i][8] != "":
                        item2 = assetslist.index(data[i][8])
                        P2 = [0] * m
                        P2[item2] = 1
                        valid = True
                    elif (
                        data[i][6] == "Classes"
                        and data[i][7] != ""
                        and data[i][8] != ""
                    ):
                        P2 = np.where(
                            asset_classes[data[i][7]].values == data[i][8], 1, 0
                        )
                        P2 = P2 / np.sum(P2)
                        valid = True
                    elif data[i][6] == "" and data[i][7] == "" and data[i][8] == "":
                        P2 = [0] * m
                        valid = True
                    if valid == True:
                        P1 = ((np.array(P1) - np.array(P2)) * d).tolist()
                        P.append(P1)
                        Q.append([data[i][5] * d])

    P = np.array(P, ndmin=2, dtype=float)
    Q = np.array(Q, ndmin=2, dtype=float)

    for i in range(len(Q)):
        if Q[i, 0] < 0:
            P[i, :] = -1.0 * P[i, :]
            Q[i, :] = -1.0 * Q[i, :]

    return P, Q


def factors_views(views, loadings, const=True):
    r"""
    Create the factors constraints matrixes C and D of the constraint
    :math:`Cw \geq D`.

    Parameters
    ----------
    constraints : DataFrame of shape (n_constraints, n_fields)
        Constraints matrix, where n_constraints is the number of constraints
        and n_fields is the number of fields of constraints matrix, the fields
        are:

        - Disabled: (bool) indicates if the constraint is enable.
        - Factor: (str) the name of the factor of the constraint.
        - Sign: (str) can be '>=' or '<='.
        - Value: (scalar) is the maximum or minimum value of the factor.

    loadings : DataFrame of shape (n_assets, n_features)
        The loadings matrix.

    Returns
    -------
    P : nd-array
        The matrix P that shows the relation among factors in each factor view.

    Q : nd-array
        The matrix Q that shows the expected return of each factor view.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------
    ::

        loadings = {'const': [0.0004, 0.0002, 0.0000, 0.0006, 0.0001, 0.0003, -0.0003],
                    'MTUM': [0.1916, 1.0061, 0.8695, 1.9996, 0.0000, 0.0000, 0.0000],
                    'QUAL': [0.0000, 2.0129, 1.4301, 0.0000, 0.0000, 0.0000, 0.0000],
                    'SIZE': [0.0000, 0.0000, 0.0000, 0.4717, 0.0000, -0.1857, 0.0000],
                    'USMV': [-0.7838, -1.6439, -1.0176, -1.4407, 0.0055, 0.5781, 0.0000],
                    'VLUE': [1.4772, -0.7590, -0.4090, 0.0000, -0.0054, -0.4844, 0.9435]}

        loadings = pd.DataFrame(loadings)

        factorsviews = {'Disabled': [False, False, False],
                        'Factor': ['MTUM', 'USMV', 'VLUE'],
                        'Sign': ['<=', '<=', '>='],
                        'Value': [0.9, -1.2, 0.3],
                        'Relative Factor': ['USMV', '', '']}

        factorsviews = pd.DataFrame(factorsviews)


    The constraint looks like this:

    .. image:: images/factorsviews.png

    It is easier to construct the constraints in excel and then upload to a
    dataframe.

    To create the matrixes P and Q we use the following command:

    ::

        P, Q = rp.factors_views(factorsviews,
                                loadings,
                                const=True)


    The matrixes P and Q looks like this:

    .. image:: images/P_fxQ_f.png

    """

    if not isinstance(views, pd.DataFrame) and not isinstance(loadings, pd.DataFrame):
        raise ValueError("constraints and loadings must be DataFrames")

    if views.shape[1] != 5:
        raise ValueError("constraints must have five columns")

    n = len(views)
    data = views.fillna("")
    data = data.values.tolist()
    factorslist = loadings.columns.tolist()
    if const == True:
        factorslist = factorslist[1:]
    m = len(factorslist)

    P = []
    Q = []
    for i in range(0, n):
        if data[i][0] == False:
            item = factorslist.index(data[i][1])
            if data[i][2] == ">=":
                d = 1
            elif data[i][2] == "<=":
                d = -1
            P1 = [0] * m
            P1[item] = d
            if data[i][4] != "":
                item = factorslist.index(data[i][4])
                P1[item] = -d
            P.append(P1)
            Q.append([data[i][3] * d])

    P = np.array(P, ndmin=2, dtype=float)
    Q = np.array(Q, ndmin=2, dtype=float)

    return P, Q


def assets_clusters(
    returns,
    custom_cov=None,
    codependence="pearson",
    linkage="ward",
    k=None,
    max_k=10,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
    leaf_order=True,
):
    r"""
    Create asset classes based on hierarchical clustering.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    linkage : string, optional
        Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
        The default is 'ward'. Possible values are:

        - 'single'.
        - 'complete'.
        - 'average'.
        - 'weighted'.
        - 'centroid'.
        - 'median'.
        - 'ward'.
        - 'DBHT'. Direct Bubble Hierarchical Tree.

    k : int, optional
        Number of clusters. This value is took instead of the optimal number
        of clusters calculated with the two difference gap statistic.
        The default is None.
    max_k : int, optional
        Max number of clusters used by the two difference gap statistic
        to find the optimal number of clusters. The default is 10.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value choice by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.
    leaf_order : bool, optional
        Indicates if the cluster are ordered so that the distance between
        successive leaves is minimal. The default is True.

    Returns
    -------
    clusters : DataFrame
        A dataframe with asset classes based on hierarchical clustering.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------

    ::

        clusters = rp.assets_clusters(returns,
                                      codependence='pearson',
                                      linkage='ward',
                                      k=None,
                                      max_k=10,
                                      alpha_tail=0.05,
                                      leaf_order=True)


    The clusters dataframe looks like this:

    .. image:: images/clusters_df.png

    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    # Calculating codependence matrix and distance metric
    codep, dist = af.codep_dist(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
    )
    # Hierarchical clustering
    dist = dist.to_numpy()
    dist = pd.DataFrame(dist, columns=codep.columns, index=codep.index)
    if linkage == "DBHT":
        # different choices for D, S give different outputs!
        D = dist.to_numpy()  # dissimilarity matrix
        if codependence in {"pearson", "spearman"}:
            S = (1 - dist**2).to_numpy()
        else:
            S = codep.copy().to_numpy()  # similarity matrix
        (_, _, _, _, _, clustering) = db.DBHTs(
            D, S, leaf_order=leaf_order
        )  # DBHT clustering
    else:
        p_dist = squareform(dist, checks=False)
        clustering = hr.linkage(p_dist, method=linkage, optimal_ordering=leaf_order)

    # Optimal number of clusters
    if k is None:
        k = af.two_diff_gap_stat(dist, clustering, max_k)

    # Building clusters
    clusters_inds = hr.fcluster(clustering, k, criterion="maxclust")
    labels = np.array(returns.columns.tolist())

    clusters = {"Assets": [], "Clusters": []}

    for i, v in enumerate(clusters_inds):
        clusters["Assets"].append(labels[i])
        clusters["Clusters"].append("Cluster " + str(v))

    clusters = pd.DataFrame(clusters)
    clusters = clusters.sort_values(by=["Assets"])

    return clusters


def hrp_constraints(constraints, asset_classes):
    r"""
    Create the upper and lower bounds constraints for hierarchical risk parity
    model.

    Parameters
    ----------
    constraints : DataFrame of shape (n_constraints, n_fields)
        Constraints matrix, where n_constraints is the number of constraints
        and n_fields is the number of fields of constraints matrix, the fields
        are:

        - Disabled: (bool) indicates if the constraint is enable.
        - Type: (str) can be: 'Assets', All Assets' and 'Each asset in a class'.
        - Position: (str) the name of the asset or asset class of the constraint.
        - Sign: (str) can be '>=' or '<='.
        - Weight: (scalar) is the maximum or minimum weight of the absolute constraint.

    asset_classes : DataFrame of shape (n_assets, n_cols)
        Asset's classes matrix, where n_assets is the number of assets and
        n_cols is the number of columns of the matrix where the first column
        is the asset list and the next columns are the different asset's
        classes sets.

    Returns
    -------
    w_max : pd.Series
        The upper bound of hierarchical risk parity weights constraints.

    w_min : pd.Series
        The lower bound of hierarchical risk parity weights constraints.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------
    ::

        asset_classes = {'Assets': ['FB', 'GOOGL', 'NTFX', 'BAC', 'WFC', 'TLT', 'SHV'],
                         'Class 1': ['Equity', 'Equity', 'Equity', 'Equity', 'Equity',
                                     'Fixed Income', 'Fixed Income'],
                         'Class 2': ['Technology', 'Technology', 'Technology',
                                     'Financial', 'Financial', 'Treasury', 'Treasury'],}

        asset_classes = pd.DataFrame(asset_classes)
        asset_classes = asset_classes.sort_values(by=['Assets'])

        constraints = {'Disabled': [False, False, False, False, False, False],
                       'Type': ['Assets', 'Assets', 'All Assets', 'All Assets',
                                'Each asset in a class', 'Each asset in a class'],
                       'Set': ['', '', '', '','Class 1', 'Class 2'],
                       'Position': ['BAC', 'FB', '', '', 'Equity', 'Treasury'],
                       'Sign': ['>=', '<=', '<=', '>=', '<=', '<='],
                       'Weight': [0.02, 0.085, 0.09, 0.01, 0.07, 0.06]}

        constraints = pd.DataFrame(constraints)

    The constraint looks like this:

    .. image:: images/HRPConstraints.png

    It is easier to construct the constraints in excel and then upload to a
    dataframe.

    To create the pd.Series w_max and w_min we use the following command:

    ::

        w_max, w_min = rp.hrp_constraints(constraints, asset_classes)


    The pd.Series w_max and w_min looks like this (all constraints were
    merged to a single upper bound for each asset):

    .. image:: images/HRP_Bounds.png

    """

    if not isinstance(constraints, pd.DataFrame) and not isinstance(
        asset_classes, pd.DataFrame
    ):
        raise ValueError("constraints and asset_classes must be DataFrames")

    if constraints.shape[1] != 6:
        raise ValueError("constraints must have six columns")

    n = len(constraints)
    data = constraints.fillna("").copy()
    assetslist = asset_classes.iloc[:, 0].values.tolist()

    w_max = pd.Series(1.0, index=assetslist)
    w_min = pd.Series(0.0, index=assetslist)

    for i in range(0, n):
        if data.loc[i, "Disabled"] == False:
            if data.loc[i, "Type"] == "Assets":
                assets = data.loc[i, "Position"]
                if data.loc[i, "Sign"] == ">=":
                    if w_min.loc[assets] <= data.loc[i, "Weight"]:
                        w_min.loc[assets] = data.loc[i, "Weight"]
                elif data.loc[i, "Sign"] == "<=":
                    if w_max.loc[assets] >= data.loc[i, "Weight"]:
                        w_max.loc[assets] = data.loc[i, "Weight"]
            elif data.loc[i, "Type"] == "All Assets":
                if data.loc[i, "Sign"] == ">=":
                    if w_min[w_min <= data.loc[i, "Weight"]].shape[0] != 0:
                        w_min[w_min <= data.loc[i, "Weight"]] = data.loc[i, "Weight"]
                elif data.loc[i, "Sign"] == "<=":
                    if w_max[w_max >= data.loc[i, "Weight"]].shape[0] != 0:
                        w_max[w_max >= data.loc[i, "Weight"]] = data.loc[i, "Weight"]
            elif data.loc[i, "Type"] == "Each asset in a class":
                label_0 = asset_classes.columns.tolist()[0]
                label_1 = data.loc[i, "Set"]
                label_2 = data.loc[i, "Position"]
                assets = asset_classes[[label_0, label_1]][
                    asset_classes[label_1] == label_2
                ]
                assets = assets["Assets"].tolist()
                if data.loc[i, "Sign"] == ">=":
                    if (
                        w_min.loc[assets][
                            w_min.loc[assets] <= data.loc[i, "Weight"]
                        ].shape[0]
                        != 0
                    ):
                        w_min.loc[assets] = np.where(
                            w_min.loc[assets] <= data.loc[i, "Weight"],
                            data.loc[i, "Weight"],
                            w_min.loc[assets],
                        )
                elif data.loc[i, "Sign"] == "<=":
                    if (
                        w_max.loc[assets][
                            w_max.loc[assets] >= data.loc[i, "Weight"]
                        ].shape[0]
                        != 0
                    ):
                        w_max.loc[assets] = np.where(
                            w_max.loc[assets] >= data.loc[i, "Weight"],
                            data.loc[i, "Weight"],
                            w_max.loc[assets],
                        )

    return w_max, w_min


def risk_constraint(asset_classes, kind="vanilla", classes_col=None):
    r"""
    Create the risk contribution constraint vector for the risk parity model.

    Parameters
    ----------
    asset_classes : DataFrame of shape (n_assets, n_cols)
        Asset's classes matrix, where n_assets is the number of assets and
        n_cols is the number of columns of the matrix where the first column
        is the asset list and the next columns are the different asset's
        classes sets. It is only used when kind value is 'classes'. The default
        value is None.

    kind : str
        Kind of risk contribution constraint vector. The default value is 'vanilla'.
        Possible values are:

        - 'vanilla': vector of equal risk contribution per asset.
        - 'classes': vector of equal risk contribution per class.

    classes_col : str or int
        If value is str, it is the column name of the set of classes from
        asset_classes dataframe. If value is int, it is the column number of
        the set of classes from asset_classes dataframe. The default
        value is None.

    Returns
    -------
    rb : nd-array
        The risk contribution constraint vector.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------
    ::

        asset_classes = {'Assets': ['FB', 'GOOGL', 'NTFX', 'BAC', 'WFC', 'TLT', 'SHV'],
                         'Class 1': ['Equity', 'Equity', 'Equity', 'Equity', 'Equity',
                                      'Fixed Income', 'Fixed Income'],
                         'Class 2': ['Technology', 'Technology', 'Technology',
                                      'Financial', 'Financial', 'Treasury', 'Treasury'],}

        asset_classes = pd.DataFrame(asset_classes)
        asset_classes = asset_classes.sort_values(by=['Assets'])
        asset_classes.reset_index(inplace=True, drop=True)

        rb = rp.risk_constraint(asset_classes
                                kind='classes',
                                classes_col='Class 1')


    """
    if not isinstance(asset_classes, pd.DataFrame):
        raise ValueError("asset_classes must be a DataFrame")

    if kind == "vanilla":
        if asset_classes.shape[1] < 1:
            raise ValueError("asset_classes must have at least one column")

        assetslist = asset_classes.iloc[:, 0].values.tolist()
        rb = np.ones((len(assetslist), 1))
        rb /= len(assetslist)

    elif kind == "classes":
        if asset_classes.shape[1] < 2:
            raise ValueError("asset_classes must have at least two columns")

        classes = asset_classes.columns.tolist()

        if isinstance(classes_col, str) and classes_col in classes:
            A = asset_classes.loc[:, classes_col].to_frame()
            col = A.columns.to_list()[0]
        elif isinstance(classes_col, int) and classes[classes_col] in classes:
            A = asset_classes.iloc[:, classes_col].to_frame()
            col = A.columns.to_list()[0]
        else:
            raise ValueError(
                "classes_col must be a valid column or column position of asset_classes"
            )

        A["rb"] = 1
        B = A.groupby([col]).count()
        A = pd.merge(A, B, left_on=col, right_index=True, how="left")
        A["rb"] = A["rb_x"] / A["rb_y"]
        A["rb"] /= A["rb"].sum()

        rb = A["rb"].to_numpy().reshape(-1, 1)

    else:
        raise ValueError(
            "The only available values for kind parameter are 'vanilla' and 'classes'"
        )

    return rb


def connection_matrix(
    returns,
    custom_cov=None,
    codependence="pearson",
    graph="MST",
    walk_size=1,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
):
    r"""
    Create a connection matrix of walks of a specific size  based on :cite:`e-Cajas10` formula..

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    graph : string, optional
        Graph used to build the adjacency matrix. The default is 'MST'.
        Possible values are:

        - 'MST': Minimum Spanning Tree.
        - 'TMFG': Plannar Maximally Filtered Graph.

    walk_size : int, optional
        Size of the walk represented by the adjacency matrix. The default is 1.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value choice by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.

    Returns
    -------
    A_p : DataFrame
        Adjacency matrix of walks of size lower and equal than 'walk_size'.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------

    ::

        A_p = rp.connection_matrix(returns,
                                   codependence="pearson",
                                   graph="MST",
                                   walk_size=1)

    The connection matrix dataframe looks like this:

    .. image:: images/Connection_df.png

    """

    if not isinstance(returns, pd.DataFrame):
        raise ValueError("returns must be a DataFrame")

    assets = returns.columns.tolist()

    # Calculating codependence matrix and distance metric
    codep, dist = af.codep_dist(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
    )

    # Adjacency Matrix Construction
    dist = dist.to_numpy()
    dist = pd.DataFrame(dist, columns=codep.columns, index=codep.index)
    if graph == "TMFG":
        # different choices for D, S give different outputs!
        D = dist.to_numpy()  # dissimilarity matrix
        if codependence in {"pearson", "spearman"}:
            S = (1 - dist**2).to_numpy()
        else:
            S = codep.copy().to_numpy()
        (_, Rpm, _, _, _, clustering) = db.DBHTs(D, S)  # DBHT clustering
        MAdj = pd.DataFrame(Rpm, index=assets, columns=assets)
        G = nx.from_pandas_adjacency(MAdj)
    elif graph == "MST":
        MAdj = nx.from_pandas_adjacency(dist)
        G = nx.minimum_spanning_tree(MAdj)
    else:
        raise ValueError("Only TMFG or MST graphs are available")

    A = nx.adjacency_matrix(G).toarray()
    A = np.where(A != 0, 1, 0)

    A_p = np.zeros_like(A)
    for i in range(walk_size + 1):
        A_p += np.linalg.matrix_power(A, i)

    n, n = A.shape
    A_p = np.clip(A_p, 0, 1) - np.identity(n)
    A_p = np.ceil(A_p)

    return A_p


def centrality_vector(
    returns,
    measure="Degree",
    custom_cov=None,
    codependence="pearson",
    graph="MST",
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
):
    r"""
    Create a centrality vector from the adjacency matrix of an asset network based on :cite:`e-Cajas10` formula.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    measure : str, optional
        Centrality measure. The default is 'Degree'. Possible values are:

        - 'Degre': Node's degree centrality. Number of edges connected to a node.
        - 'Eigenvector': Eigenvector centrality. See more in `eigenvector_centrality_numpy <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.eigenvector_centrality_numpy.html#eigenvector-centrality-numpy>`_.
        - 'Katz': Katz centrality. See more in `katz_centrality_numpy <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.katz_centrality_numpy.html#katz-centrality-numpy>`_.
        - 'Closeness': Closeness centrality. See more in `closeness_centrality <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.closeness_centrality.html#closeness-centrality>`_.
        - 'Betweeness': Betweeness centrality. See more in `betweenness_centrality <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.betweenness_centrality.html#betweenness-centrality>`_.
        - 'Communicability': Communicability betweeness centrality. See more in `communicability_betweenness_centrality <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.communicability_betweenness_centrality.html#communicability-betweenness-centrality>`_.
        - 'Subgraph': Subgraph centrality. See more in `subgraph_centrality <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.subgraph_centrality.html#subgraph-centrality>`_.
        - 'Laplacian': Laplacian centrality. See more in `laplacian_centrality <https://networkx.org/documentation/stable/reference/algorithms/generated/networkx.algorithms.centrality.laplacian_centrality.html#laplacian-centrality>`_.

    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    graph : string, optional
        Graph used to build the adjacency matrix. The default is 'MST'.
        Possible values are:

        - 'MST': Minimum Spanning Tree.
        - 'TMFG': Plannar Maximally Filtered Graph.

    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value choice by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.

    Returns
    -------
    A_p : DataFrame
        Adjacency matrix of walks of size 'walk_size'.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------

    ::

        C_v = rp.centrality_vector(returns,
                                   measure='Degree',
                                   codependence="pearson",
                                   graph="MST")

    The neighborhood matrix looks like this:

    .. image:: images/Centrality_df.png

    """

    Adj = connection_matrix(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        graph=graph,
        walk_size=1,
        bins_info=bins_info,
        gs_threshold=gs_threshold,
    )

    n, n = Adj.shape
    G = nx.from_numpy_array(Adj)
    if measure == "Degree":
        CM = np.ones((1, n)) @ Adj
    elif measure == "Eigenvector":
        CM = nx.eigenvector_centrality_numpy(G)
    elif measure == "Katz":
        CM = nx.katz_centrality_numpy(G)
    elif measure == "Closeness":
        CM = nx.closeness_centrality(G)
    elif measure == "Betweeness":
        CM = nx.betweenness_centrality(G)
    elif measure == "Communicability":
        CM = nx.communicability_betweenness_centrality(G)
    elif measure == "Subgraph":
        CM = nx.subgraph_centrality(G)
    elif measure == "Laplacian":
        CM = nx.laplacian_centrality(G)

    if measure != "Degree":
        CM = pd.Series(CM).to_numpy().reshape(1, -1)

    return CM


def clusters_matrix(
    returns,
    custom_cov=None,
    codependence="pearson",
    linkage="ward",
    k=None,
    max_k=10,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
    leaf_order=True,
):
    r"""
    Creates an adjacency matrix that represents the clusters from the hierarchical
    clustering process based on :cite:`e-Cajas11` formula.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    linkage : string, optional
        Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
        The default is 'ward'. Possible values are:

        - 'single'.
        - 'complete'.
        - 'average'.
        - 'weighted'.
        - 'centroid'.
        - 'median'.
        - 'ward'.
        - 'DBHT'. Direct Bubble Hierarchical Tree.

    k : int, optional
        Number of clusters. This value is took instead of the optimal number
        of clusters calculated with the two difference gap statistic.
        The default is None.
    max_k : int, optional
        Max number of clusters used by the two difference gap statistic
        to find the optimal number of clusters. The default is 10.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value choice by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.
    leaf_order : bool, optional
        Indicates if the cluster are ordered so that the distance between
        successive leaves is minimal. The default is True.

    Returns
    -------
    A_c : ndarray
        Adjacency matrix of clusters.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------

    ::

        C_M = rp.clusters_matrix(returns,
                                 codependence='pearson',
                                 linkage='ward',
                                 k=None,
                                 max_k=10)


    The clusters matrix looks like this:

    .. image:: images/Clusters_matrix_df.png

    """

    assets = returns.columns.tolist()
    n = len(assets)
    clusters = assets_clusters(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        linkage=linkage,
        k=k,
        max_k=max_k,
        bins_info="KN",
        alpha_tail=alpha_tail,
        gs_threshold=0.5,
        leaf_order=leaf_order,
    )

    df = pd.DataFrame([], index=assets)

    for i in clusters["Clusters"].unique():
        labels = clusters[clusters["Clusters"] == i]["Assets"].tolist()
        df1 = pd.Series(np.zeros((n,)), index=assets)
        df1[labels] = 1
        df = pd.concat([df, df1], axis=1)

    A_c = df.to_numpy()
    A_c = A_c @ A_c.T - np.identity(n)

    return A_c


def connected_assets(
    returns,
    w,
    custom_cov=None,
    codependence="pearson",
    graph="MST",
    walk_size=1,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
):
    r"""
    Calculates the percentage invested in connected assets based on :cite:`e-Cajas10` formula.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    w : DataFrame of shape (n_assets, 1)
        Portfolio weights.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    graph : string, optional
        Graph used to build the adjacency matrix. The default is 'MST'.
        Possible values are:

        - 'MST': Minimum Spanning Tree.
        - 'TMFG': Plannar Maximally Filtered Graph.

    walk_size : int, optional
        Size of the walk represented by the adjacency matrix. The default is 1.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value choice by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.

    Returns
    -------
    CA : float
        Percentage invested in connected assets.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------

    ::

        ca = rp.connected_assets(returns,
                                 w,
                                 codependence="pearson",
                                 graph="MST",
                                 walk_size=1)

    """

    A_p = connection_matrix(
        returns=returns,
        custom_cov=custom_cov,
        codependence=codependence,
        graph=graph,
        walk_size=walk_size,
        bins_info=bins_info,
        gs_threshold=gs_threshold,
    )

    n, n = A_p.shape
    ones = np.ones((n, 1))
    w_ = np.array(w)
    wwt = np.abs(w_ @ w_.T)
    ca = ones.T @ (A_p * wwt) @ ones
    ca /= ones.T @ wwt @ ones

    return ca.item()


def related_assets(
    returns,
    w,
    custom_cov=None,
    codependence="pearson",
    linkage="ward",
    k=None,
    max_k=10,
    bins_info="KN",
    alpha_tail=0.05,
    gs_threshold=0.5,
    leaf_order=True,
):
    r"""
    Calculates the percentage invested in related assets based on :cite:`e-Cajas11` formula.

    Parameters
    ----------
    returns : DataFrame
        Assets returns.
    w : DataFrame of shape (n_assets, 1)
        Portfolio weights.
    custom_cov : DataFrame or None, optional
        Custom covariance matrix, used when codependence parameter has value
        'custom_cov'. The default is None.
    codependence : str, can be {'pearson', 'spearman', 'abs_pearson', 'abs_spearman', 'distance', 'mutual_info', 'tail' or 'custom_cov'}
        The codependence or similarity matrix used to build the distance
        metric and clusters. The default is 'pearson'. Possible values are:

        - 'pearson': pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.
        - 'spearman': spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{spearman}_{i,j})}`.
        - 'kendall': kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{kendall}_{i,j})}`.
        - 'gerber1': Gerber statistic 1 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber1}_{i,j})}`.
        - 'gerber2': Gerber statistic 2 correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{gerber2}_{i,j})}`.
        - 'abs_pearson': absolute value pearson correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_spearman': absolute value spearman correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho_{i,j}|)}`.
        - 'abs_kendall': absolute value kendall correlation matrix. Distance formula: :math:`D_{i,j} = \sqrt{(1-|\rho^{kendall}_{i,j}|)}`.
        - 'distance': distance correlation matrix. Distance formula :math:`D_{i,j} = \sqrt{(1-\rho^{distance}_{i,j})}`.
        - 'mutual_info': mutual information matrix. Distance used is variation information matrix.
        - 'tail': lower tail dependence index matrix. Dissimilarity formula :math:`D_{i,j} = -\log{\lambda_{i,j}}`.
        - 'custom_cov': use custom correlation matrix based on the custom_cov parameter. Distance formula: :math:`D_{i,j} = \sqrt{0.5(1-\rho^{pearson}_{i,j})}`.

    linkage : string, optional
        Linkage method of hierarchical clustering, see `linkage <https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html?highlight=linkage#scipy.cluster.hierarchy.linkage>`_ for more details.
        The default is 'ward'. Possible values are:

        - 'single'.
        - 'complete'.
        - 'average'.
        - 'weighted'.
        - 'centroid'.
        - 'median'.
        - 'ward'.
        - 'DBHT'. Direct Bubble Hierarchical Tree.

    k : int, optional
        Number of clusters. This value is took instead of the optimal number
        of clusters calculated with the two difference gap statistic.
        The default is None.
    max_k : int, optional
        Max number of clusters used by the two difference gap statistic
        to find the optimal number of clusters. The default is 10.
    bins_info: int or str
        Number of bins used to calculate variation of information. The default
        value is 'KN'. Possible values are:

        - 'KN': Knuth's choice method. See more in `knuth_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.knuth_bin_width.html>`_.
        - 'FD': Freedman–Diaconis' choice method. See more in `freedman_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.freedman_bin_width.html>`_.
        - 'SC': Scotts' choice method. See more in `scott_bin_width <https://docs.astropy.org/en/stable/api/astropy.stats.scott_bin_width.html>`_.
        - 'HGR': Hacine-Gharbi and Ravier' choice method.
        - int: integer value choice by user.

    alpha_tail : float, optional
        Significance level for lower tail dependence index. The default is 0.05.
    gs_threshold : float, optional
        Gerber statistic threshold. The default is 0.5.
    leaf_order : bool, optional
        Indicates if the cluster are ordered so that the distance between
        successive leaves is minimal. The default is True.

    Returns
    -------
    RA : float
        Percentage invested in related assets.

    Raises
    ------
        ValueError when the value cannot be calculated.

    Examples
    --------

    ::

        ra = rp.related_assets(returns,
                               w,
                               codependence="pearson",
                               linkage="ward",
                               k=None,
                               max_k=10)

    """

    L_a = clusters_matrix(
        returns,
        custom_cov=custom_cov,
        codependence=codependence,
        linkage=linkage,
        k=k,
        max_k=max_k,
        bins_info=bins_info,
        alpha_tail=alpha_tail,
        gs_threshold=gs_threshold,
        leaf_order=leaf_order,
    )

    n, n = L_a.shape
    ones = np.ones((n, 1))
    w_ = np.array(w)
    wwt = np.abs(w_ @ w_.T)
    ra = ones.T @ (L_a * wwt) @ ones
    ra /= ones.T @ wwt @ ones

    return ra.item()
