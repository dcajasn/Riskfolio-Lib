import numpy as np
import pandas as pd


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

        import riskfolio.ConstraintsFunctions as cf

        asset_classes = {'Assets': ['FB', 'GOOGL', 'NTFX', 'BAC', 'WFC', 'TLT', 'SHV'],
                         'Class 1': ['Equity', 'Equity', 'Equity', 'Equity', 'Equity',
                                     'Fixed Income', 'Fixed Income'],
                         'Class 2': ['Technology', 'Technology', 'Technology',
                                     'Financial', 'Financial', 'Treasury', 'Treasury'],}

        asset_classes = pd.DataFrame(asset_classes)

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

        A, B = cf.assets_constraints(constraints, asset_classes)


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

    A = np.array(A, ndmin=2)
    B = np.array(B, ndmin=2)

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
                       'Value': [0.9, -1.2, 0.3],}

        constraints = pd.DataFrame(constraints)


    The constraint looks like this:

    .. image:: images/Constraints2.png

    It is easier to construct the constraints in excel and then upload to a
    dataframe.

    To create the matrixes C and D we use the following command:

    ::

        C, D = cf.factors_constraints(constraints, loadings)


    The matrixes C and D looks like this (all constraints were converted to a linear
    constraint):

    .. image:: images/CxD.png

    """

    if not isinstance(constraints, pd.DataFrame) and not isinstance(
        loadings, pd.DataFrame
    ):
        raise ValueError("constraints and loadings must be DataFrames")

    if constraints.shape[1] != 4:
        raise ValueError("constraints must have four columns")

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
            C.append(C1 * d)
            D.append([data[i][3] * d])

    C = np.array(C, ndmin=2)
    D = np.array(D, ndmin=2)

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

        P, Q = cf.assets_views(views, asset_classes)


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

    P = np.array(P, ndmin=2)
    Q = np.array(Q, ndmin=2)

    for i in range(len(Q)):
        if Q[i, 0] < 0:
            P[i, :] = -1 * P[i, :]
            Q[i, :] = -1 * Q[i, :]

    return P, Q
