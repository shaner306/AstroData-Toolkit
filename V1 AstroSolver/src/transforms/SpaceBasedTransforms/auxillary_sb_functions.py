from collections import namedtuple

from astropy.table import Table


def update_sb_final_transform_columns(sb_final_transform_columns,
                                      index,
                                      filter_fci,
                                      filter_fci_sigma,
                                      zprime_fci,
                                      zprime_fci_sigma):
    """
    Update columns to be used for the transform table based on information
    from the current image.

    Parameters
    ----------
    sb_final_transform_columns : namedtuple
        Attributes:
            index : array-like
                Name of the colour index used to calculate filter_fci and
                zprime_fci.
            filter_fci : np.float64
                T coefficient for index.
            filter_fci_sigma : np.float64
                Standard deviation of the T coefficient for index.
            zprime_fci : np.float64
                Zero point for index.
            zprime_fci_sigma : np.float64
                Standard deviation of the zero point for index.
    index : string
        Name of the colour index used to calculate filter_fci and zprime_fci.
    filter_fci : float
        T coefficient for index.
    filter_fci_sigma : float
        Standard deviation of the T coefficient for index.
    zprime_fci : float
        Zero point for index.
    zprime_fci_sigma : float
        Standard deviation of the zero point for index.

    Returns
    -------
    updated_sb_final_transform_columns : namedtuple
        Attributes:
            index : array-like
                Name of the colour index used to calculate filter_fci and
                zprime_fci.
            filter_fci : np.float64
                T coefficient for index.
            filter_fci_sigma : np.float64
                Standard deviation of the T coefficient for index.
            zprime_fci : np.float64
                Zero point for index.
            zprime_fci_sigma : np.float64
                Standard deviation of the zero point for index.

    """
    updated_sb_final_transform_columns = sb_final_transform_columns
    updated_sb_final_transform_columns.index.append(index)
    updated_sb_final_transform_columns.filter_fci.append(filter_fci)
    updated_sb_final_transform_columns.filter_fci_sigma.append(
        filter_fci_sigma)
    updated_sb_final_transform_columns.zprime_fci.append(zprime_fci)
    updated_sb_final_transform_columns.zprime_fci_sigma.append(
        zprime_fci_sigma)
    return updated_sb_final_transform_columns


def init_sb_final_transform_columns():
    """
    Initialize the columns that will create the table used for the space-based
    transforms.

    Returns
    -------
    sb_final_transform_columns : namedtuple
        Attributes:
            index : empty list
                Name of the colour index used to calculate filter_fci and
                zprime_fci.
            filter_fci : empty list
                T coefficient for index.
            filter_fci_sigma : empty list
                Standard deviation of the T coefficient for index.
            zprime_fci : empty list
                Zero point for index.
            zprime_fci_sigma : empty list
                Standard deviation of the zero point for index.
    """
    index = []
    filter_fci = []
    filter_fci_sigma = []
    zprime_fci = []
    zprime_fci_sigma = []
    sb_final_transform_columns = namedtuple('sb_final_transform_columns',
                                            ['index',
                                             'filter_fci',
                                             'filter_fci_sigma',
                                             'zprime_fci',
                                             'zprime_fci_sigma'])
    return sb_final_transform_columns(index,
                                      filter_fci,
                                      filter_fci_sigma,
                                      zprime_fci,
                                      zprime_fci_sigma)


def create_sb_final_transform_table(sb_final_transform_columns):
    """
    Convert the columns of the space-based transform table into an AstroPy
    table.

    Parameters
    ----------
    sb_final_transform_columns : namedtuple
        Attributes:
            index : array-like
                Name of the colour index used to calculate filter_fci and
                zprime_fci.
            filter_fci : np.float64
                T coefficient for index.
            filter_fci_sigma : np.float64
                Standard deviation of the T coefficient for index.
            zprime_fci : np.float64
                Zero point for index.
            zprime_fci_sigma : np.float64
                Standard deviation of the zero point for index.

    Returns
    -------
    sb_final_transform_table : astropy.table.Table
        Table containing results of the final space-based transforms.
        Has columns:
            CI : string
                Name of the colour index used to calculate the corresponding
                T_fCI and Z_fCI.
            T_fCI : float
                T coefficient for the corresponding CI.
            T_fCI_sigma : float
                Standard deviation of the T coefficient for the corresponding
                CI.
            Z_fCI : float
                Zero point for the corresponding CI.
            Z_fCI_sigma : float
                Standard deviation of the Zero point for the corresponding
                CI.

    """
    sb_final_transform_table = Table(
        names=[
            'CI',
            'T_fCI',
            'T_fCI_sigma',
            'Z_fCI',
            'Z_fCI_sigma'
        ],
        data=[
            sb_final_transform_columns.index,
            sb_final_transform_columns.filter_fci,
            sb_final_transform_columns.filter_fci_sigma,
            sb_final_transform_columns.zprime_fci,
            sb_final_transform_columns.zprime_fci_sigma
        ]
    )
    return sb_final_transform_table